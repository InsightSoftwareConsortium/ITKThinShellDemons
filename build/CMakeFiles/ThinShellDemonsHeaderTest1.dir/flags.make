# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# compile CXX with /usr/bin/c++
CXX_FLAGS =   -msse2 -O3 -DNDEBUG -fPIE  

CXX_DEFINES = 

CXX_INCLUDES = -I/home/pranjal.sahu/ITK/ITK-build/Modules/ThirdParty/Eigen3/src -I/home/pranjal.sahu/ITK/ITK-build/Modules/ThirdParty/KWIML/src -I/home/pranjal.sahu/ITK/Modules/ThirdParty/KWIML/src -I/home/pranjal.sahu/ITK/Modules/ThirdParty/VNL/src/vxl/core -I/home/pranjal.sahu/ITK/ITK-build/Modules/ThirdParty/VNL/src/vxl/v3p/netlib -I/home/pranjal.sahu/ITK/ITK-build/Modules/ThirdParty/VNL/src/vxl/vcl -I/home/pranjal.sahu/ITK/ITK-build/Modules/ThirdParty/VNL/src/vxl/core -I/home/pranjal.sahu/ITK/ITK-build/Modules/Core/Common -I/home/pranjal.sahu/ITK/Modules/Core/Common/include -I/home/pranjal.sahu/ITK/Modules/Filtering/ImageFilterBase/include -I/home/pranjal.sahu/ITK/ITK-build/Modules/ThirdParty/Netlib -I/home/pranjal.sahu/ITK/Modules/Numerics/Statistics/include -I/home/pranjal.sahu/ITK/Modules/Core/Transform/include -I/home/pranjal.sahu/ITK/Modules/Core/Mesh/include -I/home/pranjal.sahu/ITK/Modules/Core/ImageAdaptors/include -I/home/pranjal.sahu/ITK/Modules/Core/ImageFunction/include -I/home/pranjal.sahu/ITK/Modules/Filtering/ImageGrid/include -I/home/pranjal.sahu/ITK/Modules/Filtering/ImageCompose/include -I/home/pranjal.sahu/ITK/ITK-build/Modules/ThirdParty/ZLIB/src -I/home/pranjal.sahu/ITK/Modules/ThirdParty/ZLIB/src -I/home/pranjal.sahu/ITK/ITK-build/Modules/ThirdParty/MetaIO/src/MetaIO/src -I/home/pranjal.sahu/ITK/Modules/ThirdParty/MetaIO/src/MetaIO/src -I/home/pranjal.sahu/ITK/Modules/Core/SpatialObjects/include -I/home/pranjal.sahu/ITK/Modules/Filtering/ImageStatistics/include -I/home/pranjal.sahu/ITK/Modules/Filtering/Path/include -I/home/pranjal.sahu/ITK/Modules/Filtering/ImageIntensity/include -I/home/pranjal.sahu/ITK/Modules/Filtering/Smoothing/include -I/home/pranjal.sahu/ITK/Modules/Filtering/DisplacementField/include -I/home/pranjal.sahu/ITK/Modules/Numerics/Optimizers/include -I/home/pranjal.sahu/ITK/Modules/Numerics/Optimizersv4/include -I/home/pranjal.sahu/ITK/Modules/Core/FiniteDifference/include -I/home/pranjal.sahu/ITK/Modules/Filtering/ImageGradient/include -I/home/pranjal.sahu/ITK/Modules/Filtering/ImageSources/include -I/home/pranjal.sahu/ITK/Modules/Filtering/ImageFeature/include -I/home/pranjal.sahu/ITK/Modules/Registration/Common/include -I/home/pranjal.sahu/ITK/Modules/Registration/Metricsv4/include -I/home/pranjal.sahu/ITK/Modules/Core/QuadEdgeMesh/include -I/home/pranjal.sahu/ITK/Modules/Filtering/QuadEdgeMeshFiltering/include -I/home/pranjal.sahu/ITK/Modules/Registration/RegistrationMethodsv4/include -I/home/pranjal.sahu/ITK/Modules/Bridge/VTK/include -I/home/pranjal.sahu/ITK/Modules/Bridge/VtkGlue/include -I/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/include -I/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/include -isystem /home/pranjal.sahu/ITK/Modules/ThirdParty/Eigen3/src/itkeigen/.. -isystem /home/pranjal.sahu/ITK/ITK-build/Modules/ThirdParty/KWSys/src -isystem /home/pranjal.sahu/ITK/Modules/ThirdParty/VNL/src/vxl/v3p/netlib -isystem /home/pranjal.sahu/ITK/Modules/ThirdParty/VNL/src/vxl/vcl -isystem /home/pranjal.sahu/ITK/Modules/ThirdParty/VNL/src/vxl/core/vnl/algo -isystem /home/pranjal.sahu/ITK/Modules/ThirdParty/VNL/src/vxl/core/vnl -isystem /home/pranjal.sahu/ITK/ITK-build/Modules/ThirdParty/VNL/src/vxl/core/vnl -isystem /home/pranjal.sahu/VTK-9.0.3/build/IO/Image -isystem /home/pranjal.sahu/VTK-9.0.3/IO/Image -isystem /home/pranjal.sahu/VTK-9.0.3/build/Common/Core -isystem /home/pranjal.sahu/VTK-9.0.3/Common/Core -isystem /home/pranjal.sahu/VTK-9.0.3/build/Utilities/KWIML -isystem /home/pranjal.sahu/VTK-9.0.3/Utilities/KWIML -isystem /home/pranjal.sahu/VTK-9.0.3/build/Utilities/KWSys -isystem /home/pranjal.sahu/VTK-9.0.3/Utilities/KWSys -isystem /home/pranjal.sahu/VTK-9.0.3/build/Common/ExecutionModel -isystem /home/pranjal.sahu/VTK-9.0.3/Common/ExecutionModel -isystem /home/pranjal.sahu/VTK-9.0.3/build/Common/DataModel -isystem /home/pranjal.sahu/VTK-9.0.3/Common/DataModel -isystem /home/pranjal.sahu/VTK-9.0.3/build/Common/Math -isystem /home/pranjal.sahu/VTK-9.0.3/Common/Math -isystem /home/pranjal.sahu/VTK-9.0.3/build/Common/Transforms -isystem /home/pranjal.sahu/VTK-9.0.3/Common/Transforms -isystem /home/pranjal.sahu/VTK-9.0.3/build/Imaging/Core -isystem /home/pranjal.sahu/VTK-9.0.3/Imaging/Core -isystem /home/pranjal.sahu/VTK-9.0.3/build/Imaging/Sources -isystem /home/pranjal.sahu/VTK-9.0.3/Imaging/Sources -isystem /home/pranjal.sahu/VTK-9.0.3/build/Rendering/OpenGL2 -isystem /home/pranjal.sahu/VTK-9.0.3/Rendering/OpenGL2 -isystem /home/pranjal.sahu/VTK-9.0.3/build/Rendering/Core -isystem /home/pranjal.sahu/VTK-9.0.3/Rendering/Core -isystem /home/pranjal.sahu/VTK-9.0.3/build/Filters/Core -isystem /home/pranjal.sahu/VTK-9.0.3/Filters/Core -isystem /home/pranjal.sahu/VTK-9.0.3/build/Common/Misc -isystem /home/pranjal.sahu/VTK-9.0.3/Common/Misc -isystem /home/pranjal.sahu/VTK-9.0.3/build/Rendering/UI -isystem /home/pranjal.sahu/VTK-9.0.3/Rendering/UI -isystem /home/pranjal.sahu/VTK-9.0.3/build/ThirdParty/glew/vtkglew -isystem /home/pranjal.sahu/VTK-9.0.3/ThirdParty/glew/vtkglew -isystem /home/pranjal.sahu/VTK-9.0.3/build/ThirdParty/glew -isystem /home/pranjal.sahu/VTK-9.0.3/ThirdParty/glew -isystem /home/pranjal.sahu/VTK-9.0.3/build/Rendering/FreeType -isystem /home/pranjal.sahu/VTK-9.0.3/Rendering/FreeType -isystem /home/pranjal.sahu/VTK-9.0.3/build/ThirdParty/freetype/vtkfreetype -isystem /home/pranjal.sahu/VTK-9.0.3/ThirdParty/freetype/vtkfreetype -isystem /home/pranjal.sahu/VTK-9.0.3/build/ThirdParty/freetype -isystem /home/pranjal.sahu/VTK-9.0.3/ThirdParty/freetype -isystem /home/pranjal.sahu/VTK-9.0.3/build/ThirdParty/zlib/vtkzlib -isystem /home/pranjal.sahu/VTK-9.0.3/ThirdParty/zlib/vtkzlib -isystem /home/pranjal.sahu/VTK-9.0.3/build/ThirdParty/zlib -isystem /home/pranjal.sahu/VTK-9.0.3/ThirdParty/zlib -isystem /home/pranjal.sahu/VTK-9.0.3/build/Interaction/Style -isystem /home/pranjal.sahu/VTK-9.0.3/Interaction/Style -isystem /home/pranjal.sahu/VTK-9.0.3/build/Interaction/Widgets -isystem /home/pranjal.sahu/VTK-9.0.3/Interaction/Widgets -isystem /home/pranjal.sahu/VTK-9.0.3/build/Filters/General -isystem /home/pranjal.sahu/VTK-9.0.3/Filters/General -isystem /home/pranjal.sahu/VTK-9.0.3/build/Filters/Sources -isystem /home/pranjal.sahu/VTK-9.0.3/Filters/Sources 
