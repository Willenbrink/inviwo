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
    , propStepSize("stepSize", "Step size", 0.01, 0, 1)
    , propUseContrastEnhancement("useContrastEnhancement", "Use Contrast Enhancement", false)
    , propUseFastLic("useFastLic", "Use FastLIC", false)
    , propDesiredMean("DesiredMean", "Desired Mean", 255 / 2, 0, 255)
    , propDesiredStandardDeviation("DesiredStandardDeviation", "Desired Standard Deviation", 25, 0,
                                   255)
    , propKernelSize("kernelSize", "Kernel Size", 5, 1, 200, 1)
    , propColor("color", "Color", false)

{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
    // TODO: Register additional properties
    addProperty(propKernelSize);
    
    
    addProperty(propStepSize);
    addProperty(propUseContrastEnhancement);
    addProperty(propDesiredMean);
    addProperty(propDesiredStandardDeviation);
    addProperty(propUseFastLic);
    addProperty(propColor);
    
    

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
        for (int i = 0; i < propKernelSize; i++) {
            // TODO use normalized RK4. Works but which approach should we use?
            dvec2 newPoint = Integrator::RK4(vectorField, currentPoint, propStepSize, 1);
            if (!vectorField.isInside(newPoint)) {
                break;
            }
            currentPoint = newPoint;

            samples.push_back(sampleNoise(currentPoint));
        }
        currentPoint = startPoint;
        for (int i = 0; i < propKernelSize; i++) {
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
        for (int i = 0; i < propKernelSize; i++) {
            // TODO use normalized RK4. Works but which approach should we use?
            dvec2 newPoint = Integrator::RK4(vectorField, currentPoint, propStepSize, 1);
            if (!vectorField.isInside(newPoint)) {
                break;
            }
            currentPoint = newPoint;

            samples.push_back(sampleNoise(currentPoint));
        }
        currentPoint = startPoint;
        for (int i = 0; i < propKernelSize; i++) {
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

    
    int numOfNonBlackPixels = 0;
    //calculate mean
    int sum = 0;
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            if (! licImage.readPixel(size2_t(i, j)).x == 0) {
                numOfNonBlackPixels++;
                sum = sum + licImage.readPixel(size2_t(i, j)).x;
            }
            
        }
    }
    double mean = sum / numOfNonBlackPixels;

    //calculate standard deviation
    double totalDif = 0;
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            if (!licImage.readPixel(size2_t(i, j)).x == 0) {
                totalDif = totalDif + (licImage.readPixel(size2_t(i, j)).x - mean) *
                                          (licImage.readPixel(size2_t(i, j)).x - mean);
            }
        }
    }
    double variance = totalDif / numOfNonBlackPixels;
    double standardDeviation = sqrt(variance);

    
   // LogProcessorWarn("mean: " << mean); //mean: 126
    
    //LogProcessorWarn("standard devaiation: " << standardDeviation); //standard devaiation : 16.5052

    //stretching factor
    double stretchingFactor = propDesiredStandardDeviation / standardDeviation;

    //re"color"
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            if (!licImage.readPixel(size2_t(i, j)).x == 0) {
                int currentPixel = licImage.readPixel(size2_t(i, j)).x - mean;
                int newVal = propDesiredMean + stretchingFactor * (currentPixel - mean);
                licImage.setPixel(size2_t(i, j), dvec4(newVal, newVal, newVal, 255));
            }
        }
    }

    //add colors depending on velocity of vector
    if (propColor) {     
        dvec2 scale = dvec2(diff.x / texDims_.x, diff.y / texDims_.y);
        for (size_t j = 0; j < texDims_.y; j++) {                                  
            for (size_t i = 0; i < texDims_.x; i++) {
                dvec2 point = BBoxMin_ + dvec2(i * scale.x, j * scale.y);
                dvec2 vec = vectorField.interpolate(point);
                double val = glm::length(vec) * 175;
                int r, g, b;
                if (val < 128) { //it will be on the scale from red to green
                    r = 255 - 2 * val;
                    g = 2 * val;
                    b = 0;
                } else { //it will be on the scale from green to blue
                    val = val - 128;
                    r = 0;
                    g = 255 - 2 * val;
                    b = 2 * val;
                }
                int original = licImage.readPixel(size2_t(i, j)).x *0.4;

                licImage.setPixel(dvec2(i, j), dvec4(0.5 * b + original, 0.5 * g + original, 0.5 * r + original, 255)); //blue and red are swapped to look like the lecture video
                //licImage.setPixel(dvec2(i, j), dvec4(b, g, r, 255));
                //licImage.setPixel(dvec2(i, j), dvec4(std::min(b, original), std::min(g, original), std::min(r, original), 255));
            }
        }
    }
   



    // TODO contrast enhancement?

    licOut_.setData(outImage);
}

} // namespace inviwo
