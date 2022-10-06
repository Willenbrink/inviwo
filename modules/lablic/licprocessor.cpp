/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/volume/volumeram.h>
#include <lablic/licprocessor.h>
#include <labstreamlines/integrator.h>

#include <random>
#include <math.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor", // Class identifier
    "LICProcessor",            // Display name
    "KTH Labs",                // Category
    CodeState::Experimental,   // Code state
    Tags::None,                // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
    : Processor()
    , volumeIn_("volIn")
    , noiseTexIn_("noiseTexIn")
    , licOut_("licOut")
    , propMaxSteps("maxSteps", "Max Steps", 5, 0, 1000)
    , propStepSize("stepSize", "Step size", 0.01, 0, 1)
    , propUseContrastEnhancement("useContrastEnhancement", "Use Contrast Enhancement", false)
    , propUseFastLic("useFastLic", "Use FastLIC", false)
    , propDesiredMean("DesiredMean", "Desired Mean", 255 / 2, 0, 255)
    , propDesiredStandardDeviation("DesiredStandardDeviation", "Desired Standard Deviation", 25, 0, 255)
, propRandomSeed("randomSeed", "Random seed", 500, 0, 1000)


// TODO: Register additional properties
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
    // TODO: Register additional properties

    addProperty(propMaxSteps);
    addProperty(propStepSize);
    addProperty(propUseContrastEnhancement);
    addProperty(propDesiredMean);
    addProperty(propDesiredStandardDeviation);
    addProperty(propUseFastLic);
    addProperty(propRandomSeed);
    

    util::hide(propDesiredMean);
    util::hide(propDesiredStandardDeviation);

    propUseContrastEnhancement.onChange([this]() {
        if (propUseContrastEnhancement) {
            util::show(propDesiredMean);
            util::show(propDesiredStandardDeviation);
        } else {
            util::hide(propDesiredMean);
            util::hide(propDesiredStandardDeviation);
        }
    });

}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

    auto vol = volumeIn_.getData();
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    vectorFieldDims_ = vol->getDimensions();
    auto BBoxMin_ = vectorField.getBBoxMin();
    auto BBoxMax_ = vectorField.getBBoxMax();
    auto diff = BBoxMax_ - BBoxMin_;

    auto tex = noiseTexIn_.getData();
    const RGBAImage texture = RGBAImage::createFromImage(tex);
    texDims_ = tex->getDimensions();
    // Scale of tex to vectorfield
    auto scale = dvec2(diff.x / texDims_.x, diff.y / texDims_.y);

    double value = texture.readPixelGrayScale(size2_t(0, 0));

    LogProcessorInfo(value);

    // Prepare the output, it has the same dimensions as the texture and rgba values in [0,255]
    auto outImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage licImage(outImage);

    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y, 0.0));

    // Hint: Output an image showing which pixels you have visited for debugging
    std::vector<std::vector<int>> visited(texDims_.x, std::vector<int>(texDims_.y, 0));

    // TODO: Implement LIC and FastLIC
    // This code instead sets all pixels to the same gray value
    auto sampleNoise = [&](dvec2 point) {
        point -= BBoxMin_;
        point = dvec2(point.x / scale.x, point.y / scale.y);
        return texture.sample(point);
    };

    auto LIC = [&](dvec2 startPoint) {
        std::vector<dvec4> samples;
        dvec2 currentPoint = startPoint;
        samples.push_back(sampleNoise(currentPoint));
        for (int i = 0; i < propMaxSteps; i++) {
            // TODO use normalized RK4. Works but which approach should we use?
            dvec2 newPoint = Integrator::RK4(vectorField, currentPoint, propStepSize, 1);
            if (!vectorField.isInside(newPoint)) {
                break;
            }
            currentPoint = newPoint;

            samples.push_back(sampleNoise(currentPoint));
        }
        currentPoint = startPoint;
        for (int i = 0; i < propMaxSteps; i++) {
            dvec2 newPoint = Integrator::RK4(vectorField, currentPoint, propStepSize, 0);
            if (!vectorField.isInside(newPoint)) {
                break;
            }
            currentPoint = newPoint;

            samples.push_back(sampleNoise(currentPoint));
        }
        return samples;
    };

    std::map<dvec2, std::vector<dvec2>> linesMap;
    auto fastLIC = [&](dvec2 startPoint) {
        std::vector<dvec4> samples;
        dvec2 currentPoint = startPoint;
        samples.push_back(sampleNoise(currentPoint));
        for (int i = 0; i < propMaxSteps; i++) {
            // TODO use normalized RK4. Works but which approach should we use?
            dvec2 newPoint = Integrator::RK4(vectorField, currentPoint, propStepSize, 1);
            if (!vectorField.isInside(newPoint)) {
                break;
            }
            currentPoint = newPoint;

            samples.push_back(sampleNoise(currentPoint));
        }
        currentPoint = startPoint;
        for (int i = 0; i < propMaxSteps; i++) {
            dvec2 newPoint = Integrator::RK4(vectorField, currentPoint, propStepSize, 0);
            if (!vectorField.isInside(newPoint)) {
                break;
            }
            currentPoint = newPoint;

            samples.push_back(sampleNoise(currentPoint));
        }
        return samples;
    };

    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            dvec2 point = BBoxMin_ + dvec2(i * scale.x, j * scale.y);
            //TODO user-defined kernel-size
            auto samples = LIC(point);
            int val = 0;
            int len = samples.size();
            for (int c = 0; c < len; c++) {
                val += samples[c].x;
            }
            if (len)
                val = val / len;

            // int val = int(texture.readPixelGrayScale(size2_t(i, j)));
            licImage.setPixel(size2_t(i, j), dvec4(val, val, val, 255));
            // or
            // licImage.setPixelGrayScale(size2_t(i, j), val);
        }
    }

    // TODO contrast enhancement?

    licOut_.setData(outImage);
}

} // namespace inviwo
