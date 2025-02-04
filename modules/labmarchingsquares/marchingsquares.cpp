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
    , propGauss("gauss", "Gaussian Filtering")
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

    addProperty(propGauss);

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
            util::hide(propIsoValue, propIsoColor);
            util::show(propNumContours, propIsoTransferFunc);
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
    // ivec2 ij = {0, 0};
    // double valueAt00 = grid.getValueAtVertex(ij);
    LogProcessorInfo("The max coordinates are: " << bBoxMax << ".");

    // Initialize the output: mesh and vertices for the grid and bounding box
    auto gridmesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> gridvertices;

    auto indexBufferBBox = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

    // Set the random seed to the one selected in the interface
    randGenerator.seed(static_cast<std::mt19937::result_type>(propRandomSeed.get()));
    auto rand = [&]() { return randomValue(0.0, 1.0);};

    // Properties are accessed with propertyName.get()
    if (propShowGrid.get()) {
        auto indexBufferGrid = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

        // The function drawLineSegments creates two vertices at the specified positions,
        // that are placed into the Vertex vector defining our mesh.
        // An index buffer specifies which of those vertices should be grouped into to make up
        // lines/trianges/quads. Here two vertices make up a line segment.
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
    }

    // Set the created grid mesh as output
    gridmesh->addVertices(gridvertices);
    meshGridOut.setData(gridmesh);

    // Our input is const (i.e. cannot be altered), but you need to compute smoothed data and write
    // it somewhere
    ScalarField2 smoothedField = ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin);
    int sizeX = grid.getNumVerticesPerDim().x,sizeY = grid.getNumVerticesPerDim().y;
    for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
            float sum = 0.0f;
            int count = 0;
            for(int i = -2; i <= 2; i++) {
                for(int j = -2; j <= 2; j++) {
                    int multipliers[3][3] = {{15,12,5},{12,9,4},{5,4,2}};
                    int mult = multipliers[std::abs(i)][std::abs(j)];
                    if(x+i >= 0 && y+j >= 0 && x+i < sizeX && y+j < sizeY) {
                        sum += grid.getValueAtVertex({x+i,y+j}) * mult;
                        count += mult;
                    }
                }
            }

            if(propGauss.get()) {
                smoothedField.setValueAtVertex({x, y}, sum / count);
            } else {
                smoothedField.setValueAtVertex({x, y}, grid.getValueAtVertex({x,y}));
            }
        }
    }
    // Initialize the output: mesh and vertices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    auto drawIsoContour = [&](double propertyIsoValue, vec4 color) {
        for (int x = 0; x < smoothedField.getNumVerticesPerDim().x - 1; x++) {
            for (int y = 0; y < smoothedField.getNumVerticesPerDim().y - 1; y++) {
                // LogProcessorWarn("Loop through x=" << x << ",y=" << y);
                vec2 bottomLeft = vec2(x, y);
                vec2 bottomRight = vec2(x + 1, y);
                vec2 topLeft = vec2(x, y + 1);
                vec2 topRight = vec2(x + 1, y + 1);
                double bottomLeftValue = smoothedField.getValueAtVertex(bottomLeft);
                double bottomRightValue = smoothedField.getValueAtVertex(bottomRight);
                double topLeftValue = smoothedField.getValueAtVertex(topLeft);
                double topRightValue = smoothedField.getValueAtVertex(topRight);

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

                    vec2 newPos = vec2(smoothedField.getPositionAtVertex(bottomLeft).x + relative * cellSize.x,
                                       smoothedField.getPositionAtVertex(bottomLeft).y);
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

                    vec2 newPos = vec2(smoothedField.getPositionAtVertex(topLeft).x + relative * cellSize.x,
                                       smoothedField.getPositionAtVertex(topLeft).y);
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

                    vec2 newPos = vec2(smoothedField.getPositionAtVertex(bottomLeft).x,
                                       smoothedField.getPositionAtVertex(
                                           bottomLeft).y + relative * cellSize.y);
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

                    vec2 newPos = vec2(smoothedField.getPositionAtVertex(bottomRight).x,
                                       smoothedField.getPositionAtVertex(bottomRight).y + relative *
                                       cellSize.y);
                    specialPoints.push_back(newPos);

                }

                std::sort(specialPoints.begin(), specialPoints.end(), [](vec2 a, vec2 b) {
                    return a.x < b.x;
                });
                if (specialPoints.size() == 4 && rand() > 0.5f &&
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

    if (propMultiple.get() == 0) {
        double propertyIsoValue = propIsoValue.get();
        drawIsoContour(propertyIsoValue, propIsoColor.get());
    } else {
        int numContours = propNumContours.get();

        // The transfer function normalizes the input data and sampling colors
        // from the transfer function assumes normalized input, i.e. 0.0 to 1.0
        auto transfer = [&](int i) {
            return propIsoTransferFunc.get().sample(static_cast<float>(i) / numContours);
        };

        double increase = (maxValue - minValue) / numContours;
        double currentValue = minValue + increase/2;
        for (int i = 0; i< numContours; i++) {
            drawIsoContour(currentValue, transfer(i));
            currentValue = currentValue + increase;
        }
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
