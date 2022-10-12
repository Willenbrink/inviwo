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

namespace inviwo {

const vec4 Topology::ColorsCP[6] = {
    vec4(1, 1, 0, 1),    // Saddle - Yellow
    vec4(1, 0, 0, 1),    // AttractingNode - Red
    vec4(0, 0, 1, 1),    // RepellingNode - Blue
    vec4(0.5, 0, 1, 1),  // AttractingFocus - Purple
    vec4(1, 0.5, 0, 1),  // RepellingFocus - Orange
    vec4(0, 1, 0, 1)     // Center - Green
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
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
    addPort(meshBBoxOut);

    // TODO: Register additional properties
    // addProperty(propertyName);
}

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
    std::function<std::optional<dvec2>(dvec2, dvec2, int)> decomposition;
    decomposition = [&](dvec2 lowerLeft, dvec2 upperRight, int depth) -> std::optional<dvec2> {
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
        if( depth == 0 ) {
            return std::optional<dvec2>(center);
        }
        auto lowerLeft_dec = decomposition(lowerLeft, center, depth - 1);
        auto upperLeft_dec = decomposition(upperLeft, center, depth - 1);
        auto lowerRight_dec = decomposition(lowerRight, center, depth - 1);
        auto upperRight_dec = decomposition(upperRight, center, depth - 1);
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

    // Looping through all values in the vector field.
    for (size_t j = 0; j < dims[1]; ++j) {
        for (size_t i = 0; i < dims[0]; ++i) {
            dvec2 ll = to_tex(dvec2(i,j));
            dvec2 ur = to_tex(dvec2(i+1,j+1));
            auto zero = decomposition(ll, ur, 10);
            if (zero.has_value())
            {
                point(zero.value(), vec4(0,255,0,255));
            }
        }
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
    vec4 colorCenter = ColorsCP[static_cast<int>(TypeCP::Center)];

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