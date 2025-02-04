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
#include <inviwo/core/util/utilities.h>
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
    , propStepSize("stepSize", "Step size", 0.001, 0, 1)
    , propUseContrastEnhancement("useContrastEnhancement", "Use Contrast Enhancement", false)
    , propUseFastLic("useFastLic", "Use FastLIC", true)
    , propDesiredMean("DesiredMean", "Desired Mean", 0.5, 0, 1)
    , propDesiredStandardDeviation("DesiredStandardDeviation", "Desired Standard Deviation", 0.1, 0, 1)
    , propKernelSize("kernelSize", "Kernel Size", 5, 0, 200)
    , propColor("color", "Color", false)

{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
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

    // Prepare the output, it has the same dimensions as the texture and rgba values in [0,255]
    auto outImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage licImage(outImage);

    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y, 0.0));

    // This code instead sets all pixels to the same gray value
    auto vectorToPoint = [&](dvec2 point) {
        point -= BBoxMin_;
        point = dvec2(point.x / scale.x, point.y / scale.y);
        // return point;
        return vec2(point.x, point.y);
        // return vec2((int) point.x, (int) point.y);
    };
    auto sampleNoise = [&](dvec2 point) {
        return texture.sample(vectorToPoint(point)).x;
    };

    auto LIC = [&](dvec2 startPoint) {
        std::vector<double> samples;
        dvec2 currentPoint = startPoint;
        samples.push_back(sampleNoise(currentPoint));
        for (int i = 0; i < propKernelSize; i++) {
            dvec2 newPoint = Integrator::RK4_norm(vectorField, currentPoint, propStepSize, 1);
            if (!vectorField.isInside(newPoint)) {
                break;
            }
            currentPoint = newPoint;

            samples.push_back(sampleNoise(currentPoint));
        }
        currentPoint = startPoint;
        for (int i = 0; i < propKernelSize; i++) {
            dvec2 newPoint = Integrator::RK4_norm(vectorField, currentPoint, propStepSize, 0);
            if (!vectorField.isInside(newPoint)) {
                break;
            }
            currentPoint = newPoint;

            samples.push_back(sampleNoise(currentPoint));
        }
        return samples;
    };

    // Hint: Output an image showing which pixels you have visited for debugging
    // Contains the final pixel values
    std::vector<std::vector<int>> visited(texDims_.x, std::vector<int>(texDims_.y, 0));

    int maxSteps = 5000;
    auto getStreamLine = [&](dvec2 currentPoint, std::vector<dvec2> &points, vec2 texPoint, bool forward) {
        for (int steps = 0;; steps++) {
            dvec2 newPoint = Integrator::RK4_norm(vectorField, currentPoint, propStepSize, forward);
            // Stop integrating at the borders or if stuck in a loop
            if (!vectorField.isInside(newPoint) || steps >= maxSteps) {
                break;
            }
            currentPoint = newPoint;

            vec2 newTexPoint = vectorToPoint(currentPoint);
            if(newTexPoint == texPoint) {
                LogProcessorInfo("Loop detected");
                return true;
            }
            points.push_back(currentPoint);
        }
        return false;
    };

    int kernelSize = propKernelSize;
    auto fastLIC = [&](dvec2 startPoint) {
        std::vector<dvec2> points;
        vec2 texPoint = vectorToPoint(startPoint);
        if(visited[texPoint.x][texPoint.y]) {
            return;
        }
        points.push_back(startPoint);

        bool loop = getStreamLine(startPoint, points, texPoint, false);
        std::reverse(points.begin(), points.end());
        if(!loop) {
            texPoint = vectorToPoint(points[0]);
            loop = getStreamLine(startPoint, points, texPoint, true);
        }

        // Convolute
        int streamLineLenght = points.size();
        double sum = 0;
        int count = 0;
        int leftIndex = 0;
        int rightIndex = kernelSize;
        // int -> optional int to encode weird loop behaviour
        // Normally the first pixel uses only half a kernel
        // In loops it uses a full kernel
        auto getVal = [&](int i) {
            if((i >= 0 && i < streamLineLenght) || loop) {
                // Because % is not the modulo operator if negative
                return (i + streamLineLenght * maxSteps) % streamLineLenght;
            }
            return -1;
        };
        if(true) {
            //Precompute kernel for the first pixel
            // for(int i = -kernelSize - 1; i < kernelSize; i++) {
            //     int ind = getVal(i);
            //     if(ind != -1) {
            //         dvec2 point = points[ind];
            //         sum += sampleNoise(point);
            //         count++;
            //     }
            // }
            // Slide it along the points
            for(int pointIndex = 0; pointIndex < streamLineLenght; pointIndex++) {
                // int ind = getVal(pointIndex + kernelSize);
                // if(ind != -1) {
                //     dvec2 point = points[ind];
                //     sum += sampleNoise(point);
                //     count++;
                // }
                // ind = getVal(pointIndex - kernelSize - 1);
                // if(ind != -1) {
                //     dvec2 point = points[ind];
                //     sum -= sampleNoise(point);
                //     count--;
                // }
                
                if (rightIndex >= streamLineLenght - 1) {
                    sum -= sampleNoise(points[leftIndex]);
                    count--;
                    leftIndex++;
                } else if (rightIndex - leftIndex < kernelSize * 2) {
                    sum += sampleNoise(points[rightIndex]);
                    count++;
                    rightIndex++;
                } else  {
                    sum -= sampleNoise(points[leftIndex]);
                    sum += sampleNoise(points[rightIndex]);
                    leftIndex++;
                    rightIndex++;
                }
                
                vec2 point = vectorToPoint(points[pointIndex]);
                if(point.x < texDims_.x && point.y < texDims_.y && !visited[point.x][point.y])
                    visited[point.x][point.y] = sum / std::max(1,count);
            }
        } else {
            // Slide it along the points
            for(int pointIndex = 0; pointIndex < streamLineLenght; pointIndex++) {
                sum = 0;
                count = 0;
                if (leftIndex == rightIndex) {
                    break;
                }
                for(int j = leftIndex; j < rightIndex; j++) {
                    
                    // int ind = getVal(i+j);
                    // if(ind != -1) {
                    //     dvec2 point = points[ind];
                    //     sum += sampleNoise(point);
                    //     count++;
                    // }
                    dvec2 point = points[j];
                    sum += sampleNoise(point);
                    count++;
                }
                
                vec2 point = vectorToPoint(points[pointIndex]);
                if(point.x < texDims_.x && point.y < texDims_.y)
                    visited[point.x][point.y] = sum / std::max(1,count);
                
                if (rightIndex - leftIndex < kernelSize * 2) {
                    rightIndex++;
                } else if (rightIndex >= streamLineLenght - 1 ) {
                    leftIndex++;
                } else {
                    leftIndex++;
                    rightIndex++;
                }
            }
        }
    };

    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            dvec2 point = BBoxMin_ + dvec2(i * scale.x, j * scale.y);
            int val = 0;
            if(propUseFastLic) {
                fastLIC(point);
                val = visited[i][j];
            } else {
                auto samples = LIC(point);
                int len = samples.size();
                for(int c = 0; c < len; c++) {
                    val += samples[c];
                }
                if(len)
                    val = val / len;
            }

            // int val = int(texture.readPixelGrayScale(size2_t(i, j)));
            licImage.setPixel(size2_t(i, j), dvec4(val, val, val, 255));
            // or
            // licImage.setPixelGrayScale(size2_t(i, j), val);
        }
    }

    auto iterImg = [&](std::function<void(int, int, int)> f) {
        for (size_t j = 0; j < texDims_.y; j++) {
            for (size_t i = 0; i < texDims_.x; i++) {
                int val = licImage.readPixel(size2_t(i, j)).x;
                f(i,j,val);
            }
        }
    };

    if(propUseContrastEnhancement) {
        int numOfNonBlackPixels = 0;
        //calculate mean
        int sum = 0;
        iterImg([&](int, int, int val) {
            if (val) {
                sum += val;
                numOfNonBlackPixels++;
            }
        });
        double mean = sum / numOfNonBlackPixels;

        //calculate standard deviation
        double totalDif = 0;
        iterImg([&](int, int, int val) {
            if (val) {
                totalDif = totalDif + (val - mean) * (val - mean);
            }
        });
        double variance = totalDif / numOfNonBlackPixels;
        double standardDeviation = sqrt(variance);


        // LogProcessorWarn("mean: " << mean); //mean: 126

        // LogProcessorWarn("standard devaiation: " << standardDeviation); //standard devaiation : 16.5052

        //stretching factor
        double stretchingFactor = (propDesiredStandardDeviation * 256) / standardDeviation;

        //re"color"
        iterImg([&](int i, int j, int val) {
            if (val) {
                val = (propDesiredMean * 256) + stretchingFactor * (val - mean);
                licImage.setPixel(size2_t(i, j), dvec4(val, val, val, 255));
            }
        });
    }

    //add colors depending on velocity of vector
    if (propColor) {
        iterImg([&](int i, int j, int pixel) {
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
            int original = pixel * 0.4;

            licImage.setPixel(dvec2(i, j), dvec4(0.5 * b + original, 0.5 * g + original, 0.5 * r + original, 255)); //blue and red are swapped to look like the lecture video
            //licImage.setPixel(dvec2(i, j), dvec4(b, g, r, 255));
            //licImage.setPixel(dvec2(i, j), dvec4(std::min(b, original), std::min(g, original), std::min(r, original), 255));
        });
    }

    licOut_.setData(outImage);
}

} // namespace inviwo
