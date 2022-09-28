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
    , propMaxArcLenght("maxArcLenght", "Maximum arc lenght", 0.5, 0, 10)
    , propMinVelocity("minVelocity", "Minimum velocity", 0.001)
    , propUniformGrid("uniformGrid", "Use a uniform grid", true)
    , propUniformNumX("uniformNumX", "Number of sample points in X direction", 10)
    , propUniformNumY("uniformNumY", "Number of sample points in Y direction", 10)
    , propRandomMagnitude("randomMagnitude", "Prefer samples at points of high magnitude", false)
    , propNumberOfSeeds("NumberOfSeeds", "Number of random seeds", 100, 1, 1000)

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
    addProperty(propUniformGrid);
    addProperty(propUniformNumX);
    addProperty(propUniformNumY);
    addProperty(propRandomMagnitude);
    addProperty(propNumberOfSeeds);

    util::hide(propNumberOfSeeds);
    util::hide(propRandomMagnitude);
    util::hide(propUniformGrid);
    util::hide(propUniformNumX);
    util::hide(propUniformNumY);

    // Show properties for a single seed and hide properties for multiple seeds
    // (TODO)
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart, propNumStepsTaken);
            util::hide(propUniformGrid);
        } else {
            util::hide(propStartPoint, mouseMoveStart, propNumStepsTaken);
            util::show(propUniformGrid);
        }
    });
    propUniformGrid.onChange([this]() {
        if (propUniformGrid.get() == 0) {
            util::show(propNumberOfSeeds, propRandomMagnitude);
            util::hide(propUniformNumX, propUniformNumY);
        } else {
            util::show(propUniformNumX, propUniformNumY);
            util::hide(propNumberOfSeeds, propRandomMagnitude);
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

    VectorField2 smoothedField = vectorField;
    if (propNormalizeVectorField) {
        smoothedField = VectorField2(vectorField.getNumVerticesPerDim(), BBoxMin_,
                                     BBoxMax_ - BBoxMin_);
        int sizeX = vectorField.getNumVerticesPerDim().x, sizeY = vectorField.getNumVerticesPerDim()
                .y;
        for (int x = 0; x < sizeX; x++) {
            for (int y = 0; y < sizeY; y++) {
                auto val = vectorField.getValueAtVertex({x, y});
                val = glm::normalize(val);
                smoothedField.setValueAtVertex({x, y}, val);
            }
        }
    }
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferStreamLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

    auto calcStreamline = [this, &indexBufferPoints, &indexBufferStreamLines, &smoothedField, &vertices, &red](vec2 startPoint) {
        // Draw start point
        if (propDisplayPoints.get() != 0)
            Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

        // TODO: Create one stream line from the given start point
        dvec2 currentPoint = startPoint;
        double arcLength = 0;
        int i = 0;
        for (; i < propMaxSteps && arcLength < propMaxArcLenght; i++) {
            dvec2 newPoint = Integrator::RK4(smoothedField, currentPoint, propStepSize,
                                             propDirection == 0);
            dvec2 movement = newPoint - currentPoint;
            if (!smoothedField.isInside(newPoint)
                || glm::length(movement) < propMinVelocity) {
                break;
            }
            double distance = sqrt((newPoint.x - currentPoint.x) * (newPoint.x - currentPoint.x) +
                                   (newPoint.y - currentPoint.y) * (newPoint.y - currentPoint.y));
            arcLength = + distance;
            Integrator::drawLineSegment(currentPoint, newPoint, red, indexBufferStreamLines.get(),
                                        vertices);
            if (propDisplayPoints)
                Integrator::drawPoint(newPoint, red, indexBufferPoints.get(),
                                      vertices);
            currentPoint = newPoint;

            dvec2 value = smoothedField.interpolate(currentPoint);
            if (std::abs(value.x) < FLT_EPSILON && std::abs(value.y) < FLT_EPSILON) {
                break;
            }
        }

        // TODO: Use the propNumStepsTaken property to show how many steps have actually been
        // integrated This could be different from the desired number of steps due to stopping
        // conditions (too slow, boundary, ...)
        propNumStepsTaken.set(i);
    };

    if (propSeedMode.get() == 0) {
        vec2 startPoint = propStartPoint.get();
        calcStreamline(startPoint);
    } else {
        if (propUniformGrid.get() == 0) {
            for (int j = 0; j < propNumberOfSeeds; j++) {
                double drm = double(RAND_MAX);
                float x = (rand() - drm/2) / drm * (BBoxMax_.x - BBoxMin_.x);
                float y = (rand() - drm/2) / drm * (BBoxMax_.y - BBoxMin_.y);

                vec2 startPoint = vec2(x, y);
                calcStreamline(startPoint);
            }
        } else {
            // (TODO: Bonus, sample randomly according to magnitude of the vector field)
            int max_x = propUniformNumX.get();
            int max_y = propUniformNumY.get();
            const dvec2 stepSize = {(BBoxMax_.x - BBoxMin_.x) / max_x, (BBoxMax_.y - BBoxMin_.y) / max_y};
            for(int i = 0; i < max_x; i++) {
                for(int j = 0; j < max_y; j++) {
                    dvec2 startPoint = BBoxMin_ + dvec2(stepSize.x * i, stepSize.y * j);
                    calcStreamline(startPoint);
                }
            }
        }
    }

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

} // namespace inviwo
