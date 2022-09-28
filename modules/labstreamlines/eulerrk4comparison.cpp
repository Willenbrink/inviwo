/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp
 *  Init    : Tuesday, September 19, 2017 - 15:08:24
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/algorithm/boundingbox.h>
#include <inviwo/core/interaction/events/mouseevent.h>
#include <labstreamlines/eulerrk4comparison.h>
#include <labstreamlines/integrator.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo EulerRK4Comparison::processorInfo_{
    "org.inviwo.EulerRK4Comparison", // Class identifier
    "Euler RK4 Comparison",          // Display name
    "KTH Lab",                       // Category
    CodeState::Experimental,         // Code state
    Tags::None,                      // Tags
};

const ProcessorInfo EulerRK4Comparison::getProcessorInfo() const { return processorInfo_; }

EulerRK4Comparison::EulerRK4Comparison()
    : Processor()
    , inData("inData")
    , meshOut("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(0), vec2(1024), vec2(0.5))
    , mouseMoveStart(
        "mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
        MouseButton::Left, MouseState::Press | MouseState::Move)
    ,
    // TODO: Initialize additional properties
    propIntegrationStepsEuler("integrationSteps", "Integration steps Euler", 100)
    , propStepSizeEuler("stepSizeEuler", "Step size", 0.1f)
    , propIntegrationStepsRungeKutta("integrationStepsRungeKutta", "Integration steps Runge-Kutta",
                                     10)
    , propStepSizeRungeKutta("stepSizeRungeKutta", "Step size Runge-Kutta", 0.5f)
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(meshOut);
    addPort(meshBBoxOut);
    addPort(inData);

    // Register Properties
    addProperty(propStartPoint);
    addProperty(mouseMoveStart);

    // TODO: Register additional properties
    addProperty(propIntegrationStepsEuler);
    addProperty(propStepSizeEuler);
    addProperty(propIntegrationStepsRungeKutta);
    addProperty(propStepSizeRungeKutta);

}

template <typename T, typename S> 
        inline dvec2 operator*(dvec2 const & v, float const & s)
{
    return dvec2(
    v.x * s,
    v.y * s);
}

void EulerRK4Comparison::eventMoveStart(Event* event) {
    if (!inData.hasData()) return;
    auto mouseEvent = static_cast<MouseEvent*>(event);
    vec2 mousePos = mouseEvent->posNormalized();

    // Map to bounding box range
    mousePos[0] *= static_cast<float>(BBoxMax_[0] - BBoxMin_[0]);
    mousePos[1] *= static_cast<float>(BBoxMax_[1] - BBoxMin_[1]);
    mousePos += static_cast<vec2>(BBoxMin_);

    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

void EulerRK4Comparison::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();

    // The start point should be inside the volume (set maximum to the upper right corner)
    propStartPoint.setMinValue(BBoxMin_ - dvec2(1, 1));
    propStartPoint.setMaxValue(BBoxMax_ + dvec2(1, 1));

    // Initialize mesh, vertices and index buffers for the two streamlines and the points
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    auto indexBufferEuler = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;

    // Make bounding box without vertex duplication, instead of line segments which duplicate
    // vertices, create line segments between each added points with connectivity type of the index
    // buffer
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    vec4 red = vec4(1, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin_[0], BBoxMax_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax_[0], BBoxMin_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    // Draw start point
    dvec2 startPoint = propStartPoint.get();
    Integrator::drawPoint(startPoint, black, indexBufferPoints.get(), vertices);

    // TODO: Implement the Euler and Runge-Kutta of 4th order integration schemes

    //Euler 
    dvec2 currentPoint = startPoint;
    for (int i = 0; i < propIntegrationStepsEuler; i++) {
        if (!vectorField.isInside(currentPoint)) {
            break;
        }
        auto direction = vectorField.interpolate(currentPoint);
        auto newPoint = currentPoint + direction.operator*=(propStepSizeEuler);
        Integrator::drawLineSegment(currentPoint, newPoint, black, indexBufferEuler.get(),
                                    vertices);
        currentPoint = newPoint;
    }

    // Runge-Kutta

    currentPoint = startPoint;
    for (int i = 0; i < propIntegrationStepsRungeKutta; i++) {
        if (!vectorField.isInside(currentPoint)) {
            break;
        }

        dvec2 v1 = vectorField.interpolate(currentPoint);
        dvec2 v2 = vectorField.interpolate(currentPoint + dvec2(v1.x * (propStepSizeRungeKutta / 2.0f), v1.y * (propStepSizeRungeKutta / 2.0f)));
        dvec2 v3 = vectorField.interpolate(currentPoint + dvec2(v2.x * (propStepSizeRungeKutta / 2), v2.y * (propStepSizeRungeKutta / 2)));
        dvec2 v4 = vectorField.interpolate(currentPoint + dvec2(v3.x * propStepSizeRungeKutta, v3.y * propStepSizeRungeKutta));
        dvec2 huh = v1 / 6 + v2 / 3 + v3/ 3 + v4 / 6;
        dvec2 newPoint = currentPoint + dvec2(huh.x * propStepSizeRungeKutta, huh.y * propStepSizeRungeKutta);

        Integrator::drawLineSegment(currentPoint, newPoint, red, indexBufferRK.get(),
                                    vertices);
        currentPoint = newPoint;
    }

    // and then integrate forward for a specified number of integration steps and a given stepsize
    // (these should be additional properties of the processor)

    // Integrator::Euler(vectorField, startPoint, ...);
    // Integrator::Rk4(vectorField, dims, startPoint, ...);

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

} // namespace inviwo
