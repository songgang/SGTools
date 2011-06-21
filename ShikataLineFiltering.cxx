/*=========================================================================
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkHessianToShikataMeasureImageFilter.h"
#include "itkShikataMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"

#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkShiftScaleInPlaceImageFilter.h"

#include "antsCommandLineParser.h"

int itkMultiScaleHessianBasedMeasureImageFilterTest(int argc, char *argv[]);

int main(int argc, char *argv[]) {

    //    itk::MultiThreader::SetGlobalMaximumNumberOfThreads( 1 );


    return itkMultiScaleHessianBasedMeasureImageFilterTest(argc, argv);
}

template<class ParserPointerType>
void InitializeOption(ParserPointerType &parser) {
    typedef typename ParserPointerType::ObjectType ParserType;
    typedef typename ParserType::OptionType OptionType;

    {
        typename OptionType::Pointer option = OptionType::New();
        option->SetLongName(std::string("sigma-list"));
        option->SetDescription(
                std::string(
                        "sigma list for multi-scale line filtering, eg. 1x1.414x2.x2.828"));
        option->AddValue(std::string("1x1.414x2x2.828x4"));
        parser->AddOption(option);
    }

    {
        typename OptionType::Pointer option = OptionType::New();
        option->SetLongName(std::string("distance-threshold-list"));
        option->SetDescription(
                std::string(
                        "distance threshold list for multi-scale line filtering, coupled with sigma-list"
                            "eg. 1x1x2x2"));
        option->AddValue(std::string("0x1x1x2x2.828"));
        parser->AddOption(option);
    }

    {
        typename OptionType::Pointer option = OptionType::New();

        option->SetLongName(std::string("distance-map-image"));
        option->SetDescription(
                std::string(
                        "output file name of distance map image of input mask "));
        option->AddValue(std::string(""));
        parser->AddOption(option);
    }

    {
        typename OptionType::Pointer option = OptionType::New();

        option->SetLongName(std::string("hessian-image"));
        option->SetDescription(
                std::string("output file name of hessian map image,"
                    "for debug use, each pixel is the hessian from "
                    "the selected scale"));
        option->AddValue(std::string(""));
        parser->AddOption(option);
    }

    {
        typename OptionType::Pointer option = OptionType::New();

        option->SetLongName(std::string("scale-image"));
        option->SetDescription(std::string("output file name of scale image,"
            "each pixel is the select scale from "
            "sigma-list"));
        option->AddValue(std::string(""));
        parser->AddOption(option);
    }

}

template<class TVec, class TStream>
void PrintVector(TVec &v, TStream &s){
    unsigned int n = v.size();
    if (n==0) return;
    s << v[0];
    for(unsigned int i=1; i<n;i++)
        s << "x" << v[i];
    s << std::endl;
}

int itkMultiScaleHessianBasedMeasureImageFilterTest(int argc, char *argv[]) {

    typedef itk::ants::CommandLineParser ParserType;
    ParserType::Pointer parser = ParserType::New();
    InitializeOption(parser);

    if (argc < 3) {
        std::cerr << "Missing Parameters: " << argv[0]
                << "InputScalarImage InputMaskImage LineFilteringOutputImage"
                << std::endl;
        parser->PrintMenu(std::cerr, 5, false);
        //                << "[--sigma-list AxBxC]"
        //                << "[--distance-threshold-list DxExF]"
        //                << "[--distance-map-image DistanceMapImage]"
        //                << "[--hessian-image HessianImageFileName]"
        //                << "[--scale-image ScaleImageFileName]" << std::endl;

        //
        //                << " EnhancedOutputImage ScalesOutputImage "
        //                << " [SigmaMin SigmaMax NumberOfScales ObjectDimension Bright/Dark HessianImage]"
        //                << std::endl;
        return EXIT_FAILURE;
    }
    //const double SQRT2 = vnl_math::sqrt2; //sqrt(2.0);
    std::vector<double> sigmaList;
//    sigmaList.resize(5);
//    sigmaList[0] = 1;
//    sigmaList[1] = SQRT2;
//    sigmaList[2] = 2;
//    sigmaList[3] = 2 * SQRT2;
//    sigmaList[4] = 4;

    std::vector<double> distanceThresholdList;
//    distanceThresholdList.resize(5);
//    distanceThresholdList[0] = 0;
//    distanceThresholdList[1] = 1;
//    distanceThresholdList[2] = 1;
//    distanceThresholdList[3] = 2;
//    distanceThresholdList[4] = 2 * SQRT2;

    //    std::vector<double> sigmaList;
    //    sigmaList.resize(1);
    //    sigmaList[0] = 2*SQRT2;
    //
    //    std::vector<double> distanceThresholdList;
    //    distanceThresholdList.resize(1);
    //    distanceThresholdList[0] = 0;


    parser->Parse(argc - 4, argv + 4);

    std::string inputScalarImageFileName = argv[1];
    std::string inputMaskImageFileName = argv[2];
    std::string outputLineFilteredImageFileName = argv[3];

    sigmaList = parser->ConvertVector<double> (
            parser->GetOption("sigma-list")->GetValue());
    distanceThresholdList = parser->ConvertVector<double> (
            parser->GetOption("distance-threshold-list")->GetValue());

    std::cout << "sigma-list: " << std::endl;
    PrintVector(sigmaList, std::cout);

    std::cout << "distance-threshold-list: " << std::endl;
    PrintVector(distanceThresholdList, std::cout);



    std::string outputDistanceMapImageFileName = parser->GetOption(
            "distance-map-image")->GetValue();
    std::string outputHessianImageFileName =
            parser->GetOption("hessian-image")->GetValue();
    std::string outputScaleImageFileName =
            parser->GetOption("scale-image")->GetValue();

    //    char *inputScalarImageFileName = argv[1];
    //    char *inputMaskImageFileName = argv[2];
    //    char *outputLineFilteredImageFileName = argv[3];
    //    char *outputDistanceMapImageFileName = argv[4];
    //    char *outputHessianImageFileName = argv[5];
    //    char *outputScaleImageFileName = argv[6];


    // Define the dimension of the images
    const unsigned int Dimension = 3;

    typedef float InputPixelType;
    typedef itk::Image<InputPixelType, Dimension> InputImageType;

    typedef unsigned char MaskPixelType;
    typedef itk::Image<MaskPixelType, Dimension> MaskImageType;

    typedef itk::Image<float, Dimension> DistanceImageType;

    typedef float OutputPixelType;
    typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

    typedef itk::ImageFileReader<InputImageType> FileReaderType;

    typedef itk::ImageFileReader<MaskImageType> MaskFileReaderType;

    typedef itk::ImageFileWriter<OutputImageType> FileWriterType;

    typedef itk::NumericTraits<InputPixelType>::RealType RealPixelType;

    typedef itk::SymmetricSecondRankTensor<RealPixelType, Dimension>
            HessianPixelType;
    typedef itk::Image<HessianPixelType, Dimension> HessianImageType;

    // Declare the type of enhancement filter
    typedef itk::HessianToShikataMeasureImageFilter<HessianImageType,
            InputImageType, OutputImageType> ObjectnessFilterType;

    // Declare the type of multiscale enhancement filter
    typedef itk::ShikataMultiScaleHessianBasedMeasureImageFilter<
            InputImageType, HessianImageType, ObjectnessFilterType,
            OutputImageType> MultiScaleEnhancementFilterType;

    FileReaderType::Pointer imageReader = FileReaderType::New();
    imageReader->SetFileName(inputScalarImageFileName);
    try {
        imageReader->Update();
    } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
    }

    MaskFileReaderType::Pointer maskReader = MaskFileReaderType::New();
    maskReader->SetFileName(inputMaskImageFileName);
    try {
        maskReader->Update();
    } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
    }

    ObjectnessFilterType::Pointer objectnessFilter =
            ObjectnessFilterType::New();
    objectnessFilter->SetScaleObjectnessMeasure(false);
    objectnessFilter->SetBrightObject(true);

    //generate the distance image here
    typedef itk::SignedDanielssonDistanceMapImageFilter<MaskImageType,
            DistanceImageType> DistanceMapFilterType;

    DistanceMapFilterType::Pointer distanceMapFilter =
            DistanceMapFilterType::New();

    distanceMapFilter->SetInput(maskReader->GetOutput());
    distanceMapFilter->Update();

    typedef itk::ShiftScaleInPlaceImageFilter<DistanceImageType>
            ScaleFilterType;
    ScaleFilterType::Pointer scaleFilter = ScaleFilterType::New();
    scaleFilter->SetInput(distanceMapFilter->GetOutput());
    scaleFilter->SetScale(-1.0);
    scaleFilter->SetShift(0.0);
    scaleFilter->Update();

    DistanceImageType::Pointer distanceImage = scaleFilter->GetOutput();

    MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter =
            MultiScaleEnhancementFilterType::New();
    multiScaleEnhancementFilter->SetInputScalarImage(imageReader->GetOutput());

    multiScaleEnhancementFilter->SetInputDistanceImage(distanceImage);

    multiScaleEnhancementFilter->SetSigmaList(sigmaList);
    multiScaleEnhancementFilter->SetDistanceThresholdList(distanceThresholdList);

    multiScaleEnhancementFilter->SetHessianToMeasureFilter(objectnessFilter);
    //  multiScaleEnhancementFilter->SetSigmaStepMethodToLogarithmic();
    //Change NonNegativeHessianBasedMeasure to Off and regnerate vesselness image
    multiScaleEnhancementFilter->NonNegativeHessianBasedMeasureOff();

    if (outputScaleImageFileName.compare("") != 0)
        multiScaleEnhancementFilter->SetGenerateScalesOutput(true);

    if (outputHessianImageFileName.compare("") != 0)
        multiScaleEnhancementFilter->SetGenerateHessianOutput(true);

    try {
        multiScaleEnhancementFilter->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr << e << std::endl;
    }

    FileWriterType::Pointer writer = FileWriterType::New();
    writer->SetFileName(outputLineFilteredImageFileName);
    writer->UseCompressionOn();
    writer->SetInput(multiScaleEnhancementFilter->GetOutput());

    try {
        writer->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr << e << std::endl;
    }

    if (outputScaleImageFileName.compare("") != 0) {

        std::cout << "output scale image:" << outputScaleImageFileName << std::endl;
        writer->SetFileName(outputScaleImageFileName);
        writer->UseCompressionOn();
        writer->SetInput(multiScaleEnhancementFilter->GetScalesOutput());

        try {
            writer->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr << e << std::endl;
        }
    }

    if (outputDistanceMapImageFileName.compare("") != 0) {
        std::cout << "output distance map:" << outputDistanceMapImageFileName
                << std::endl;
        writer->SetFileName(outputDistanceMapImageFileName);
        writer->UseCompressionOn();
        writer->SetInput(distanceImage);

        try {
            writer->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr << e << std::endl;
        }

    }

    if (outputHessianImageFileName.compare("") != 0) {
        const HessianImageType * hessianImage =
                multiScaleEnhancementFilter->GetHessianOutput();

        std::cout << "write hessian image to "
                            << outputHessianImageFileName << std::endl;
        // write the hessian image
        typedef itk::ImageFileWriter<HessianImageType> HessianFileWriterType;
        HessianFileWriterType::Pointer writer3 = HessianFileWriterType::New();
        writer3->SetFileName(outputHessianImageFileName);
        writer3->UseCompressionOn();
        writer3->SetInput(hessianImage);

        try {

            writer3->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr << e << std::endl;
        }

    }

    return EXIT_SUCCESS;
}
