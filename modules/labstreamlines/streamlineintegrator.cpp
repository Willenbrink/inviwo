/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/interaction/events/mouseevent.h>
#include <inviwo/core/util/utilities.h>
#include <labstreamlines/integrator.h>
#include <labstreamlines/streamlineintegrator.h>
#include <labutils/scalarvectorfield.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming
// scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator", // Class identifier
    "Streamline Integrator",           // Display name
    "KTH Lab",                         // Category
    CodeState::Experimental,           // Code state
    Tags::None,                        // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const { return processorInfo_; }

StreamlineIntegrator::StreamlineIntegrator()
    : Processor()
    , inData("volIn")
    , meshOut("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , propDisplayPoints("displayPoints", "Display Points", true)
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(-1), vec2(1), vec2(0.1))
    , propSeedMode("seedMode", "Seeds")
    , propNumStepsTaken("numstepstaken", "Number of actual steps", 0, 0, 100000)
    , mouseMoveStart(
        "mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
        MouseButton::Left, MouseState::Press | MouseState::Move)
    , propDirection("direction", "Integration direction")
    , propStepSize("stepSize", "Step size", 0.5)
    , propNormalizeVectorField("normalizeVectorField", "Normalize vector field", false)
    , propMaxSteps("maxSteps", "Maximum number of steps", 18)
    , propMaxArcLenght("maxArcLenght", "Maximum arc lenght", 100, 0, 1000)
    , propMinVelocity("minVelocity", "Minimum velocity", 0.001)
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional),
// increment (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(inData);
    addPort(meshOut);
    addPort(meshBBoxOut);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("multiple", "Multiple Seeds", 1);
    addProperty(propSeedMode);
    addProperty(propStartPoint);
    addProperty(propDisplayPoints);
    addProperty(propNumStepsTaken);
    propNumStepsTaken.setReadOnly(true);
    propNumStepsTaken.setSemantics(PropertySemantics::Text);
    propDirection.addOption("left", "Forward", 0);
    propDirection.addOption("right", "Backward", 1);
    addProperty(mouseMoveStart);
    addProperty(propDirection);
    addProperty(propStepSize);
    addProperty(propNormalizeVectorField);
    addProperty(propMaxSteps);
    addProperty(propMaxArcLenght);
    addProperty(propMinVelocity);

    // Show properties for a single seed and hide properties for multiple seeds
    // (TODO)
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart, propNumStepsTaken);
            // util::hide(...)
        } else {
            util::hide(propStartPoint, mouseMoveStart, propNumStepsTaken);
            // util::show(...)
        }
    });
}

void StreamlineIntegrator::eventMoveStart(Event* event) {
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

void StreamlineIntegrator::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    auto vectorField = VectorField2::createFieldFromVolume(vol);
    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();

    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;

    // Make bounding box without vertex duplication, instead of line segments which duplicate
    // vertices, create line segments between each added points with connectivity type of the index
    // buffer
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    vec4 red = vec4(1, 0, 0, 1);
    vec4 blue = vec4(0, 0, 1, 1);
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

    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    if (propSeedMode.get() == 0) {
        auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
        auto indexBufferStreamLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
        vec2 startPoint = propStartPoint.get();
        // Draw start point
        if (propDisplayPoints.get() != 0)
            Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

        // TODO: Create one stream line from the given start point
        vec2 currentPoint = startPoint;
        int i = 0;
        for (; i < propMaxSteps; i++) {
            dvec2 newPoint = Integrator::RK4(vectorField, currentPoint, 0.5f, propDirection == 0);
            if (!vectorField.isInside(newPoint)) {
                break;
            }
            Integrator::drawLineSegment(currentPoint, newPoint, red, indexBufferStreamLines.get(),
                                        vertices);
            Integrator::drawPoint(newPoint, red, indexBufferPoints.get(), vertices);
            currentPoint = newPoint;

            dvec2 value = vectorField.interpolate(currentPoint);
            if (!vectorField.isInside(currentPoint)
                || std::abs(value.x) < 0.01 && std::abs(value.y) < 0.01) {
                break;
            }
        }

        // TODO: Use the propNumStepsTaken property to show how many steps have actually been
        // integrated This could be different from the desired number of steps due to stopping
        // conditions (too slow, boundary, ...)
        propNumStepsTaken.set(i);

    } else {
        // TODO: Seed multiple stream lines either randomly or using a uniform grid
        // (TODO: Bonus, sample randomly according to magnitude of the vector field)
    }

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

} // namespace inviwo
