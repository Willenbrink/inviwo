/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp, Anke Friederici
 *  Init    : Monday, September 11, 2017 - 12:58:42
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/util/utilities.h>
#include <math.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo MarchingSquares::processorInfo_{
    "org.inviwo.MarchingSquares", // Class identifier
    "Marching Squares",           // Display name
    "KTH Lab",                    // Category
    CodeState::Experimental,      // Code state
    Tags::None,                   // Tags
};

const ProcessorInfo MarchingSquares::getProcessorInfo() const { return processorInfo_; }


MarchingSquares::MarchingSquares()
    : Processor()
    , inData("volumeIn")
    , meshIsoOut("meshIsoOut")
    , meshGridOut("meshGridOut")
    , propShowGrid("showGrid", "Show Grid")
    , propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f), vec4(0.0f),
                    vec4(1.0f), vec4(0.1f), InvalidationLevel::InvalidOutput,
                    PropertySemantics::Color)
    , propDeciderType("deciderType", "Decider Type")
    , propRandomSeed("seed", "Random Seed", 0, 0, std::mt19937::max())
    , propMultiple("multiple", "Iso Levels")
    , propIsoValue("isovalue", "Iso Value")
    , propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f), vec4(0.0f), vec4(1.0f),
                   vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
    , propIsoTransferFunc("isoTransferFunc", "Colors", &inData) {
    // Register ports
    addPort(inData);
    addPort(meshIsoOut);
    addPort(meshGridOut);

    // Register properties
    addProperty(propShowGrid);
    addProperty(propGridColor);

    addProperty(propDeciderType);
    propDeciderType.addOption("asymptotic", "Asymptotic", 0);
    propDeciderType.addOption("random", "Random", 1);

    addProperty(propRandomSeed);
    propRandomSeed.setSemantics(PropertySemantics::Text);

    addProperty(propMultiple);

    propMultiple.addOption("single", "Single", 0);
    addProperty(propIsoValue);
    addProperty(propIsoColor);

    propMultiple.addOption("multiple", "Multiple", 1);
    addProperty(propNumContours);
    addProperty(propIsoTransferFunc);

    // The default transfer function has just two blue points
    propIsoTransferFunc.get().clear();
    propIsoTransferFunc.get().add(0.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.get().add(1.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.setCurrentStateAsDefault();

    util::hide(propGridColor, propRandomSeed, propNumContours, propIsoTransferFunc);

    propDeciderType.onChange([this]() {
        if (propDeciderType.get() == 1) {
            util::show(propRandomSeed);
        } else {
            util::hide(propRandomSeed);
        }
    });

    // Show the grid color property only if grid is actually displayed
    propShowGrid.onChange([this]() {
        if (propShowGrid.get()) {
            util::show(propGridColor);
        } else {
            util::hide(propGridColor);
        }
    });

    // Show options based on display of one or multiple iso contours
    propMultiple.onChange([this]() {
        if (propMultiple.get() == 0) {
            util::show(propIsoValue, propIsoColor);
            util::hide(propNumContours, propIsoTransferFunc);
        } else {
            util::hide(propIsoValue);
            util::show(propIsoColor, propNumContours);

            // TODO (Bonus): Comment out above if you are using the transfer function
            // and comment in below instead
            // util::hide(propIsoValue, propIsoColor);
            // util::show(propNumContours, propIsoTransferFunc);
        }
    });
}


void MarchingSquares::process() {
    if (!inData.hasData()) {
        return;
    }

    

    // Create a structured grid from the input volume
    auto vol = inData.getData();
    auto grid = ScalarField2::createFieldFromVolume(vol);

    // Extract the minimum and maximum value from the input data
    const double minValue = grid.getMinValue();
    const double maxValue = grid.getMaxValue();

    // Set the range for the isovalue to that minimum and maximum
    propIsoValue.setMinValue(minValue);
    propIsoValue.setMaxValue(maxValue);

    // You can print to the Inviwo console with Log-commands:
    LogProcessorInfo("This scalar field contains values between " << minValue << " and " << maxValue
        << ".");
    // You can also inform about errors and warnings:
    // LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
    // LogProcessorError("I am letting you know about an error"); // Will print error message in red
    // (There is also LogNetwork...() and just Log...(), these display a different source,
    // LogProcessor...() for example displays the name of the processor in the workspace while
    // Log...() displays the identifier of the processor (thus with multiple processors of the
    // same kind you would not know which one the information is coming from

    // Get the definition of our structured grid with
    // - number of vertices in each dimension {nx, ny}
    const ivec2 nVertPerDim = grid.getNumVerticesPerDim();
    // - bounding box {xmin, ymin} - {xmax, ymax}
    const dvec2 bBoxMin = grid.getBBoxMin();
    const dvec2 bBoxMax = grid.getBBoxMax();
    // - cell size {dx, dy}
    const dvec2 cellSize = grid.getCellSize();

    // Values at the vertex positions can be accessed by the indices of the vertex
    // with index i ranging between [0, nx-1] and j in [0, ny-1]
    ivec2 ij = {0, 0};
    double valueAt00 = grid.getValueAtVertex(ij);
    LogProcessorInfo("The max coordinates are: " << bBoxMax << ".");

    // Initialize the output: mesh and vertices for the grid and bounding box
    auto gridmesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> gridvertices;

    auto indexBufferBBox = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // bottomLeft to topLeft
    // drawLineSegment(bBoxMin, vec2(bBoxMin[0], bBoxMax[1]), propGridColor.get(),
    //                 indexBufferBBox.get(), gridvertices);
    // // topLeft to topRight
    // drawLineSegment(vec2(bBoxMin[0], bBoxMax[1]), bBoxMax, propGridColor.get(),
    //                 indexBufferBBox.get(), gridvertices);
    // // topRight to bottomRight
    // drawLineSegment(bBoxMax, vec2(bBoxMax[0], bBoxMin[1]), propGridColor.get(),
    //                 indexBufferBBox.get(), gridvertices);
    // // bottomRight to bottomLeft
    // drawLineSegment(vec2(bBoxMax[0], bBoxMin[1]), bBoxMin, propGridColor.get(),
    //                 indexBufferBBox.get(), gridvertices);

    // Set the random seed to the one selected in the interface
    randGenerator.seed(static_cast<std::mt19937::result_type>(propRandomSeed.get()));
    // You can create a random sample between min and max with
    float minRand = 0.0;
    float maxRand = 1.0;
    float rand = randomValue(minRand, maxRand);
    LogProcessorInfo("The first random sample for seed " << propRandomSeed.get() << " between "
        << minRand << " and " << maxRand << " is "
        << rand << ".");

    // Properties are accessed with propertyName.get()
    if (propShowGrid.get()) {
        // TODO: Add grid lines of the given color
        auto indexBufferGrid = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

        for (int x = 0; x < grid.getNumVerticesPerDim().x; x++) {

            vec2 point1 = grid.getPositionAtVertex(vec2(x, 0));
            vec2 point2 = grid.getPositionAtVertex(vec2(x, grid.getNumVerticesPerDim().y - 1));
            drawLineSegment(point1, point2, propGridColor.get(), indexBufferGrid.get(),
                            gridvertices);

        }

        for (int y = 0; y < grid.getNumVerticesPerDim().y; y++) {
            vec2 point1 = grid.getPositionAtVertex(vec2(0, y));
            vec2 point2 = grid.getPositionAtVertex(vec2(grid.getNumVerticesPerDim().x - 1, y));

            drawLineSegment(point1, point2, propGridColor.get(), indexBufferGrid.get(),
                            gridvertices);
        }
        // The function drawLineSegments creates two vertices at the specified positions,
        // that are placed into the Vertex vector defining our mesh.
        // An index buffer specifies which of those vertices should be grouped into to make up
        // lines/trianges/quads. Here two vertices make up a line segment.

        // Draw a line segment from v1 to v2 with a the given color for the grid
        // vec2 v1 = vec2(0.5, 0.5);
        // vec2 v2 = vec2(0.7, 0.7);
        // drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), gridvertices);

    }

    // Set the created grid mesh as output
    gridmesh->addVertices(gridvertices);
    meshGridOut.setData(gridmesh);

    // TODO (Bonus) Gaussian filter
    // Our input is const (i.e. cannot be altered), but you need to compute smoothed data and write
    // it somewhere
    // Create an editable structured grid with ScalarField2 smoothedField =
    // ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin); Values can be set with
    // smoothedField.setValueAtVertex({0, 0}, 4.2);
    // and read again in the same way as before
    // smoothedField.getValueAtVertex(ij);
    // Initialize the output: mesh and vertices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    LogProcessorWarn("HUH?");
    auto drawIsoContour = [&](double propertyIsoValue, double gridSize, ScalarField2 grid, vec4 color) {
        for (int x = 0; x < grid.getNumVerticesPerDim().x - 1; x++) {
            for (int y = 0; y < grid.getNumVerticesPerDim().y - 1; y++) {
                // LogProcessorWarn("Loop through x=" << x << ",y=" << y);
                vec2 bottomLeft = vec2(x, y);
                vec2 bottomRight = vec2(x + 1, y);
                vec2 topLeft = vec2(x, y + 1);
                vec2 topRight = vec2(x + 1, y + 1);
                double bottomLeftValue = grid.getValueAtVertex(bottomLeft);
                double bottomRightValue = grid.getValueAtVertex(bottomRight);
                double topLeftValue = grid.getValueAtVertex(topLeft);
                double topRightValue = grid.getValueAtVertex(topRight);

                std::vector<vec2> specialPoints;

                double diffBottomLeft = bottomLeftValue - propertyIsoValue;
                double diffBottomRight = bottomRightValue - propertyIsoValue;
                double diffTopLeft = topLeftValue - propertyIsoValue;
                double diffTopRight = topRightValue - propertyIsoValue;

                if (signbit(diffBottomLeft) != signbit(diffBottomRight)) {
                    double low, high;
                    bool highIsLeft;
                    if (bottomLeftValue < bottomRightValue) {
                        low = bottomLeftValue;
                        high = bottomRightValue;
                        highIsLeft = false;
                    } else {
                        high = bottomLeftValue;
                        low = bottomRightValue;
                        highIsLeft = true;
                    }

                    double diff = high - low;
                    double middle = propertyIsoValue - low;
                    double relative = middle / diff;
                    if (highIsLeft) {
                        relative = 1 - relative;
                    }

                    vec2 newPos = vec2(grid.getPositionAtVertex(bottomLeft).x + relative * gridSize,
                                       grid.getPositionAtVertex(bottomLeft).y);
                    specialPoints.push_back(newPos);

                }

                if (signbit(diffTopLeft) != signbit(diffTopRight)) {
                    double low, high;
                    bool highIsLeft;
                    if (topLeftValue < topRightValue) {
                        low = topLeftValue;
                        high = topRightValue;
                        highIsLeft = false;
                    } else {
                        high = topLeftValue;
                        low = topRightValue;
                        highIsLeft = true;
                    }

                    double diff = high - low;
                    double middle = propertyIsoValue - low;
                    double relative = middle / diff;
                    if (highIsLeft) {
                        relative = 1 - relative;
                    }

                    vec2 newPos = vec2(grid.getPositionAtVertex(topLeft).x + relative * gridSize,
                                       grid.getPositionAtVertex(topLeft).y);
                    specialPoints.push_back(newPos);

                }

                if (signbit(diffBottomLeft) != signbit(diffTopLeft)) {
                    double low, high;
                    bool highIsTop;
                    if (bottomLeftValue < topLeftValue) {
                        low = bottomLeftValue;
                        high = topLeftValue;
                        highIsTop = false;
                    } else {
                        high = bottomLeftValue;
                        low = topLeftValue;
                        highIsTop = true;
                    }

                    double diff = high - low;
                    double middle = propertyIsoValue - low;
                    double relative = middle / diff;
                    if (highIsTop) {
                        relative = 1 - relative;
                    }

                    vec2 newPos = vec2(grid.getPositionAtVertex(bottomLeft).x,
                                       grid.getPositionAtVertex(
                                           bottomLeft).y + relative * gridSize);
                    specialPoints.push_back(newPos);
                }

                if (signbit(diffBottomRight) != signbit(diffTopRight)) {
                    double low, high;
                    bool highIsTop;
                    if (bottomRightValue < topRightValue) {
                        low = bottomRightValue;
                        high = topRightValue;
                        highIsTop = false;
                    } else {
                        high = bottomRightValue;
                        low = topRightValue;
                        highIsTop = true;
                    }

                    double diff = high - low;
                    double middle = propertyIsoValue - low;
                    double relative = middle / diff;
                    if (highIsTop) {
                        relative = 1 - relative;
                    }

                    vec2 newPos = vec2(grid.getPositionAtVertex(bottomRight).x,
                                       grid.getPositionAtVertex(bottomRight).y + relative *
                                       gridSize);
                    specialPoints.push_back(newPos);

                }

                std::sort(specialPoints.begin(), specialPoints.end(), [](vec2 a, vec2 b) {
                    return a.x < b.x;
                });
                if (specialPoints.size() == 4 && randomValue(minRand, maxRand) > 0.5f &&
                    propDeciderType.get() == 1) {
                    auto tmp = specialPoints[1];
                    specialPoints[1] = specialPoints[2];
                    specialPoints[2] = tmp;
                }

                for (int i = 0; i < specialPoints.size(); i += 2) {
                    drawLineSegment(specialPoints[i], specialPoints[i + 1], color,
                                    indexBufferLines.get(), vertices);
                }
            }
        }
    };
    double gridSize = grid.getPositionAtVertex(vec2(1, 0)).x - grid.getPositionAtVertex(
                              vec2(0, 0)).x;

    if (propMultiple.get() == 0) {
        
        // TODO: Draw a single isoline at the specified isovalue (propIsoValue)
        // and color it with the specified color (propIsoColor)

        
        double propertyIsoValue = propIsoValue.get();

        LogProcessorWarn("Grid size: " << gridSize);

        drawIsoContour(propertyIsoValue, gridSize, grid, propIsoColor.get());

    } else {
        // TODO: Draw the given number (propNumContours) of isolines between
        // the minimum and maximum value
        int numContours = propNumContours.get();

        double increase = (maxValue - minValue) / numContours;
        double currentValue = minValue + increase/2;
        for (int i = 0; i< numContours; i++) {
            vec3 color = (vec3) (propIsoColor.get() * (numContours-i) / numContours)
                + vec3(1.0,0.0,0.0) * (i) / numContours;
            drawIsoContour(currentValue, gridSize, grid, vec4(color,1.0));
            currentValue = currentValue + increase;
        }

        // TODO (Bonus): Use the transfer function property to assign a color
        // The transfer function normalizes the input data and sampling colors
        // from the transfer function assumes normalized input, that means
        // vec4 color = propIsoTransferFunc.get().sample(0.0f);
        // is the color for the minimum value in the data
        // vec4 color = propIsoTransferFunc.get().sample(1.0f);
        // is the color for the maximum value in the data
    }

    // Note: It is possible to add multiple index buffers to the same mesh,
    // thus you could for example add one for the grid lines and one for
    // each isoline
    // Also, consider to write helper functions to avoid code duplication
    // e.g. for the computation of a single iso contour

    mesh->addVertices(vertices);
    meshIsoOut.setData(mesh);
}


float MarchingSquares::randomValue(const float min, const float max) const {
    return min + uniformReal(randGenerator) * (max - min);
}

void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
                                      IndexBufferRAM* indexBuffer,
                                      std::vector<BasicMesh::Vertex>& vertices) {
    // Add first vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    // A vertex has a position, a normal, a texture coordinate and a color
    // we do not use normal or texture coordinate, but still have to specify them
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    // Add second vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

} // namespace inviwo
