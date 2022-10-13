/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 **********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labstreamlines/integrator.h>
#include <labutils/scalarvectorfield.h>
#include <labtopo/topology.h>
#include <labtopo/utils/gradients.h>
#include <cstddef>
#include <optional>
#include "utils/gradients.h"

namespace inviwo {

const vec4 Topology::ColorsCP[7] = {
    vec4(1, 1, 0, 1),    // Saddle - Yellow
    vec4(1, 0, 0, 1),    // AttractingNode - Red
    vec4(0, 0, 1, 1),    // RepellingNode - Blue
    vec4(0.5, 0, 1, 1),  // AttractingFocus - Purple
    vec4(1, 0.5, 0, 1),  // RepellingFocus - Orange
    vec4(0, 1, 0, 1),    // Center - Green
    vec4(0, 0, 0, 1)     // Border - Black
};

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",    // Class identifier
    "Vector Field Topology",  // Display name
    "KTH Lab",                // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const { return processorInfo_; }

Topology::Topology()
    : Processor()
    , inData("inData")
    , outMesh("meshOut")
    , meshBBoxOut("meshBBoxOut")
// TODO: Initialize additional properties
    , thresholdProp("tresh", "Center Treshold", 0.0001, 0, 1)
    , stepsProp("steps", "Number of Steps", 1000, 0, 10000)
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
    addPort(meshBBoxOut);
    addProperty(thresholdProp);
    addProperty(stepsProp);

    // TODO: Register additional properties
    // addProperty(propertyName);
}

void calcBorderSwitchPoints(dvec2 start, dvec2 end, function<bool(dvec2)> isPointingInside,
                            function<bool(dvec2)> foundBorder, int maxIterations,
                            function<void(dvec2, vec4)> point, function<void(dvec2, bool)> drawStreamline) {
    if (maxIterations <= 0) {
        point(start, Topology::ColorsCP[(int) Topology::TypeCP::Border]);
        drawStreamline(start, true);
        drawStreamline(start, false);
        return;
    }
    if (isPointingInside(start) == isPointingInside(end)) {
        return;
    }
    if (foundBorder(start)) {
        // draw
        point(start, Topology::ColorsCP[(int) Topology::TypeCP::Border]);
        drawStreamline(start, true);
        drawStreamline(start, false);
        return;
    } else if (foundBorder(end)) {
        // draw
        point(end, Topology::ColorsCP[(int) Topology::TypeCP::Border]);
        drawStreamline(start, true);
        drawStreamline(start, false);
        return;
    }

    auto center = dvec2((start.x + end.x) / 2, (start.y + end.y) / 2);
    // point(center, Topology::ColorsCP[(int) Topology::TypeCP::Border]);
    
    if (foundBorder(center)) {
        //draw
        point(center, Topology::ColorsCP[(int) Topology::TypeCP::Border]);
        drawStreamline(start, true);
        drawStreamline(start, false);
        return;
    }

    if (isPointingInside(start) != isPointingInside(center)) {
        // start to center
        calcBorderSwitchPoints(start, center, isPointingInside, foundBorder, maxIterations - 1, point, drawStreamline);
    } else if (isPointingInside(center) != isPointingInside(end)) {
        calcBorderSwitchPoints(center, end, isPointingInside, foundBorder, maxIterations - 1, point, drawStreamline);
    }
};

void Topology::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);

    // Add a bounding box to the mesh
    const dvec2& BBoxMin = vectorField.getBBoxMin();
    const dvec2& BBoxMax = vectorField.getBBoxMax();
    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin[0], BBoxMax[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax[0], BBoxMin[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    size2_t dims = vectorField.getNumVerticesPerDim();
    auto diff = BBoxMax - BBoxMin;
    auto scale = dvec2(diff.x / dims.x, diff.y / dims.y);
    auto to_tex = [&](dvec2 point) {
        return dvec2(point.x * scale.x, point.y * scale.y) + BBoxMin;
    };
    auto to_vector = [&](dvec2 point) {
        point -= BBoxMin;
        return dvec2(point.x / scale.x, point.y / scale.y);
    };

    // Initialize mesh, vertices and index buffers for seperatrices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer, two consecutive points
    // make up one line), or use several index buffers with connectivity type strip.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines,
    // ConnectivityType::Strip);

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.

    auto point = [&](dvec2 v, vec4 color) {
        // drawLineSegment(start, end, color, indexBufferPoints.get(), vertices);
        indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
        vertices.push_back({dvec3(v[0], v[1], 0), dvec3(0, 0, 1), dvec3(v[0], v[1], 0), color});
    };

    // lower left, upper right, lower right, upper left
    std::function<std::optional<dvec2>(dvec2, dvec2)> decomposition;
    decomposition = [&](dvec2 lowerLeft, dvec2 upperRight) -> std::optional<dvec2> {
        dvec2 upperLeft = dvec2(lowerLeft.x, upperRight.y);
        dvec2 lowerRight = dvec2(upperRight.x, lowerLeft.y);

        dvec2 lowerLeft_v = vectorField.interpolate(lowerLeft);
        dvec2 upperLeft_v = vectorField.interpolate(upperLeft);
        dvec2 lowerRight_v = vectorField.interpolate(lowerRight);
        dvec2 upperRight_v = vectorField.interpolate(upperRight);

        int x_pos_num = 0
            + (lowerLeft_v.x > 0)
            + (upperLeft_v.x > 0)
            + (lowerRight_v.x > 0)
            + (upperRight_v.x > 0)
            ;
        int y_pos_num = 0
            + (lowerLeft_v.y > 0)
            + (upperLeft_v.y > 0)
            + (lowerRight_v.y > 0)
            + (upperRight_v.y > 0)
            ;
        if( (x_pos_num == 0 || x_pos_num == 4) || (y_pos_num == 0 || y_pos_num == 4) ) {
            return std::nullopt;
        }
        dvec2 diff = upperRight - lowerLeft;
        dvec2 center = dvec2(lowerLeft.x + diff.x / 2, lowerLeft.y + diff.y / 2);
        // point(center, vec4(0,20,0,255));
        if( std::abs(diff.x) <= 4*FLT_EPSILON && std::abs(diff.y) <= 4*FLT_EPSILON ) {
            return std::optional<dvec2>(center);
        }
        auto lowerLeft_dec = decomposition(lowerLeft, center);
        auto upperLeft_dec = decomposition(upperLeft, center);
        auto lowerRight_dec = decomposition(lowerRight, center);
        auto upperRight_dec = decomposition(upperRight, center);
        if(lowerLeft_dec.has_value()) {
            return lowerLeft_dec;
        }
        if(upperLeft_dec.has_value()) {
            return upperLeft_dec;
        }
        if(lowerRight_dec.has_value()) {
            return lowerRight_dec;
        }
        if(upperRight_dec.has_value()) {
            return upperRight_dec;
        }
        return std::nullopt;
    };

    

    auto calcStreamline = [&](vec2 startPoint, bool forward) {
        dvec2 currentPoint = startPoint;
        double arcLength = 0;
        int i = 0;
        for (; i < stepsProp; i++) {
            dvec2 newPoint = Integrator::RK4(vectorField, currentPoint, 0.01, forward);
            dvec2 movement = newPoint - currentPoint;
            if (!vectorField.isInside(newPoint)) {
                break;
            }
            double distance = sqrt((newPoint.x - currentPoint.x) * (newPoint.x - currentPoint.x) +
                                   (newPoint.y - currentPoint.y) * (newPoint.y - currentPoint.y));
            arcLength += distance;
            drawLineSegment(currentPoint, newPoint, vec4(1,1,1,1), indexBufferSeparatrices.get(), vertices);
            currentPoint = newPoint;

            dvec2 value = vectorField.interpolate(currentPoint);
            if (std::abs(value.x) < FLT_EPSILON && std::abs(value.y) < FLT_EPSILON) {
                break;
            }
        }
    };

    // Looping through all values in the vector field.
    for (size_t j = 0; j < dims[1]; ++j) {
        for (size_t i = 0; i < dims[0]; ++i) {
            dvec2 ll = to_tex(dvec2(i,j));
            dvec2 ur = to_tex(dvec2(i+1,j+1));
            auto zero = decomposition(ll, ur);
            if (zero.has_value())
            {
                dvec2 crit = zero.value();
                vec4 color = vec4(1,1,1,1);
                mat2 jacobian = vectorField.derive(crit);
                double det = jacobian[0][0] * jacobian[1][1] - jacobian[1][0] * jacobian[0][1];
                if(det == 0) {
                    continue;
                }
                util::EigenResult eigen = util::eigenAnalysis(jacobian);
                vec2 real = eigen.eigenvaluesRe;
                vec2 img = eigen.eigenvaluesIm;

                if(std::abs(img.x) <= 4*FLT_EPSILON && std::abs(img.y) <= 4*FLT_EPSILON) {
                    if(real.x > 0 && real.y > 0) {
                        color = ColorsCP[(int) TypeCP::RepellingNode];
                    }
                    if(real.x < 0 && real.y < 0) {
                        color = ColorsCP[(int) TypeCP::AttractingNode];
                    }
                    if(real.x * real.y < 0) {
                        color = ColorsCP[(int) TypeCP::Saddle];
                        mat2 eig_vec = eigen.eigenvectors;
                        dvec2 v1 = dvec2(eig_vec[0][0], eig_vec[0][1]);
                        dvec2 v2 = dvec2(eig_vec[1][0], eig_vec[1][1]);
                        calcStreamline(crit + v1 * 0.001, (real.x > 0));
                        calcStreamline(crit - v1 * 0.001, (real.x > 0));
                        calcStreamline(crit + v2 * 0.001, (real.y > 0));
                        calcStreamline(crit - v2 * 0.001, (real.y > 0));
                    }
                } else {
                    if(real.x > 0 && real.y > 0) {
                        color = ColorsCP[(int) TypeCP::RepellingFocus];
                    }
                    if(real.x < 0 && real.y < 0) {
                        color = ColorsCP[(int) TypeCP::AttractingFocus];
                    }
                    if(std::abs(real.x) <= thresholdProp && std::abs(real.y) <= thresholdProp) {
                        color = ColorsCP[(int) TypeCP::Center];
                    }
                }
                point(crit, color);
            }
        }
    }

    // Looping through the boundaries in the vector field
    

    auto isPointingLeft = [&](dvec2 location) {
        return vectorField.interpolate(location).x < 0;
    };
    auto isPointingRight = [&](dvec2 location) {
        return vectorField.interpolate(location).x > 0;
    };
    auto isPointingDown = [&](dvec2 location) {
        return vectorField.interpolate(location).y < 0;
    };
    auto isPointingUp = [&](dvec2 location) {
        return vectorField.interpolate(location).y > 0;
    };
    auto xIsZero = [&](dvec2 location) {
        return abs(vectorField.interpolate(location).x) < FLT_EPSILON;
    };
    auto yIsZero = [&](dvec2 location) {
        return abs(vectorField.interpolate(location).y) < FLT_EPSILON;
    };
    
    dvec2 bottomLeft = to_tex(dvec2(0, 0));
    dvec2 bottomRight = to_tex(dvec2(dims[0], 0));
    dvec2 topLeft = to_tex(dvec2(0, dims[1]));
    dvec2 topRight = to_tex(dvec2(dims[0], dims[1]));
    double yDiff = topLeft.y - bottomLeft.y;
    double xDiff = bottomRight.x - bottomLeft.x;
    int maxIterations = 20;
    int amountOfSteps = 10;

    // bottom left to bottom right
    double y = bottomLeft.y;
    auto lastLocation = bottomLeft;
    auto lastPointingInside = isPointingUp(lastLocation);
    for (double x = topLeft.x; x <= topRight.x; x += xDiff / amountOfSteps) {
        auto location = dvec2(x, y);
        auto value = vectorField.interpolate(location);
        auto pointingInside = isPointingUp(location);

        if (lastPointingInside != pointingInside) {
            calcBorderSwitchPoints(lastLocation, location, isPointingUp, xIsZero, maxIterations, point, calcStreamline);
        }

        lastPointingInside = pointingInside;
        lastLocation = location;
    }

    // top left to top right
    y = topLeft.y;
    lastLocation = topLeft;
    lastPointingInside = isPointingDown(lastLocation);
    for (double x = topLeft.x; x <= topRight.x; x += xDiff / amountOfSteps) {
        auto location = dvec2(x,y);
        auto value = vectorField.interpolate(location);
        auto pointingInside = isPointingDown(location);
        // LogProcessorWarn("At " << location.x << " Vector is pointing inside: " << pointingInside)
        
        if (lastPointingInside != pointingInside) {
            
            // point(location, ColorsCP[(int) Topology::TypeCP::Border] );
            calcBorderSwitchPoints(lastLocation, location, isPointingDown, xIsZero, maxIterations, point, calcStreamline);
        }
        lastPointingInside = pointingInside;
        lastLocation = location;
    }

    // bottom left to top left
    double x = bottomLeft.x;
    lastLocation = bottomLeft;
    lastPointingInside = isPointingRight(lastLocation);
    for (double y = bottomLeft.y; y <= topLeft.y; y += yDiff / amountOfSteps) {
        auto location = dvec2(x, y);
        auto value = vectorField.interpolate(location);
        auto pointingInside = isPointingRight(location);
        LogProcessorWarn("At " << location.x << ", " << location.y << " Vector is pointing inside: " << pointingInside)
        // point(location, Topology::ColorsCP[(int) Topology::TypeCP::Border]);
        
        if (lastPointingInside != pointingInside) {
            
            // point(location, ColorsCP[(int) Topology::TypeCP::Border] );
            calcBorderSwitchPoints(lastLocation, location, isPointingRight, yIsZero, maxIterations, point, calcStreamline);
        }
        lastPointingInside = pointingInside;
        lastLocation = location;
    }

    // bottom right to top right
    x = bottomLeft.x;
    lastLocation = bottomLeft;
    lastPointingInside = isPointingLeft(lastLocation);
    for (double y = bottomRight.y; y <= topRight.y; y += yDiff / amountOfSteps) {
        auto location = dvec2(x,y);
        auto value = vectorField.interpolate(location);
        auto pointingInside = isPointingLeft(location);
        // LogProcessorWarn("At " << location.x << " Vector is pointing inside: " << pointingInside)
        
        if (lastPointingInside != pointingInside) {
            
            // point(location, ColorsCP[(int) Topology::TypeCP::Border] );
            calcBorderSwitchPoints(lastLocation, location, isPointingLeft, yIsZero, maxIterations, point, calcStreamline);
        }
        lastPointingInside = pointingInside;
        lastLocation = location;
    }
    
    // Other helpful functions
    // dvec2 pos = vectorField.getPositionAtVertex(size2_t(i, j));
    // Computing the jacobian at a position
    // dmat2 jacobian = vectorField.derive(pos);
    // Doing the eigen analysis
    // auto eigenResult = util::eigenAnalysis(jacobian);
    // The result of the eigen analysis has attributed eigenvaluesRe eigenvaluesIm and
    // eigenvectors

    // Accessing the colors
    // vec4 colorCenter = ColorsCP[static_cast<int>(TypeCP::Center)];

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

void Topology::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                               IndexBufferRAM* indexBuffer,
                               std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

}  // namespace inviwo

/*

    a)
    X = y + 3
    Y = -x + 5

    b)
    X = x + 2
    Y = -(y - 7) = -y + 7


    c)
    X = (y + 3) * (x + 2)
    Y = (-x + 5) * (-y + 7)

    d)
    X = sin(y)
    Y = cos(x)

*/