# [THINSHELLDEMONS] Read Mesh and Create TSD object

# Build the Module with Python Wrapping and import the itk package which is available under Wrapping/Generators/Python
#import sys
#sys.path.append('/home/pranjal.sahu/ITKPR7/ITK/ITK-build/Wrapping/Generators/Python')

import itkConfig
itkConfig.LazyLoading = False

import numpy as np
import itk
from itk import itkThinShellDemonsMetricv4Python as tsd
from itk import itkConjugateGradientLineSearchOptimizerv4Python as itkc
import math


fixedMesh = itk.meshread('../test/Baseline/fixedMesh.vtk', itk.D)
movingMesh = itk.meshread('../test/Baseline/movingMesh.vtk', itk.D)

fixedMesh.BuildCellLinks()
movingMesh.BuildCellLinks()

PixelType = itk.D
Dimension = 3

MeshType = itk.Mesh[itk.D, Dimension]
FixedImageType = itk.Image[PixelType, Dimension]
MovingImageType = itk.Image[PixelType, Dimension]


# For getting the Bounding Box
ElementIdentifierType = itk.UL
CoordType = itk.F
Dimension = 3

VecContType = itk.VectorContainer[
    ElementIdentifierType, itk.Point[CoordType, Dimension]
]
bounding_box = itk.BoundingBox[ElementIdentifierType, Dimension, CoordType, VecContType].New()
bounding_box.SetPoints(movingMesh.GetPoints())
bounding_box.ComputeBoundingBox()

minBounds = np.array(bounding_box.GetMinimum())
maxBounds = np.array(bounding_box.GetMaximum())

imageDiagonal = 100
spacing = np.sqrt(bounding_box.GetDiagonalLength2()) / imageDiagonal
diff = maxBounds - minBounds

fixedImageSize = [0]*3
fixedImageSize[0] = math.ceil( 1.2 * diff[0] / spacing )
fixedImageSize[1] = math.ceil( 1.2 * diff[1] / spacing )
fixedImageSize[2] = math.ceil( 1.2 * diff[2] / spacing )

fixedImageOrigin = [0]*3
fixedImageOrigin[0] = minBounds[0] - 0.1 * diff[0]
fixedImageOrigin[1] = minBounds[1] - 0.1 * diff[1]
fixedImageOrigin[2] = minBounds[2] - 0.1 * diff[2]

fixedImageSpacing = np.ones(3)*spacing
fixedImageDirection = np.identity(3)


fixedImage = FixedImageType.New()
fixedImage.SetRegions(fixedImageSize)
fixedImage.SetOrigin( fixedImageOrigin )
fixedImage.SetDirection( fixedImageDirection )
fixedImage.SetSpacing( fixedImageSpacing )
fixedImage.Allocate()

VectorType = itk.Vector[itk.D, Dimension]
zeroVector = VectorType()
zeroVector.Fill(0)

DisplacementFieldType = itk.Image[VectorType, Dimension]

displacementField = DisplacementFieldType.New()
displacementField.CopyInformation(fixedImage)
displacementField.SetRegions(fixedImage.GetBufferedRegion())
displacementField.Allocate()
displacementField.FillBuffer(zeroVector)

inverseDisplacementField = DisplacementFieldType.New()
inverseDisplacementField.CopyInformation(fixedImage)
inverseDisplacementField.SetRegions(fixedImage.GetBufferedRegion())
inverseDisplacementField.Allocate()
inverseDisplacementField.FillBuffer(zeroVector)

TransformType = itk.DisplacementFieldTransform[itk.D, Dimension]

DisplacementFieldRegistrationType = itk.SyNImageRegistrationMethod[FixedImageType, MovingImageType,
            TransformType, FixedImageType, MeshType]
registration  = DisplacementFieldRegistrationType.New()

OutputTransformType = TransformType
outputTransform = OutputTransformType.New()
outputTransform.SetDisplacementField(displacementField)
outputTransform.SetInverseDisplacementField(inverseDisplacementField)
registration.SetInitialTransform(outputTransform)
registration.InPlaceOn()


# Affine Transform Object
AffineTransformType = itk.AffineTransform.D3
transform = AffineTransformType.New()
transform.SetIdentity()
transform.SetCenter(minBounds + (maxBounds - minBounds)/2)

print('Transform Created')
print(transform)


MetricType = tsd.itkThinShellDemonsMetricv4MD3
metric = MetricType.New()
metric.SetStretchWeight(1)
metric.SetBendWeight(5)
metric.SetGeometricFeatureWeight(10)
metric.UseConfidenceWeightingOn()
metric.UseMaximalDistanceConfidenceSigmaOn()
metric.UpdateFeatureMatchingAtEachIterationOff()
metric.SetMovingTransform(transform)
# Reversed due to using points instead of an image
# to keep semantics the same as in itkThinShellDemonsTest.cxx
# For the ThinShellDemonsMetricv4 the fixed mesh is regularized
metric.SetFixedPointSet(movingMesh)
metric.SetMovingPointSet(fixedMesh)
metric.SetVirtualDomainFromImage(fixedImage)
metric.Initialize()

print('TSD Metric Created')

varianceForUpdateField = spacing*spacing*25
varianceForTotalField = 0.0
registration.SetGaussianSmoothingVarianceForTheUpdateField(varianceForUpdateField)
registration.SetGaussianSmoothingVarianceForTheTotalField(varianceForTotalField)
registration.SetFixedPointSet(movingMesh)
registration.SetMovingPointSet(fixedMesh)
registration.SetMovingInitialTransform(transform)
registration.SetMetric(metric)


numberOfLevels = 3
registration.SetNumberOfLevels(numberOfLevels)
numberOfIterationsPerLevel = itk.Array[itk.UI]()
numberOfIterationsPerLevel.SetSize(numberOfLevels)
numberOfIterationsPerLevel[0] = 5
numberOfIterationsPerLevel[1] = 10
numberOfIterationsPerLevel[2] = 50
registration.SetNumberOfIterationsPerLevel(numberOfIterationsPerLevel)

shrinkFactorsPerLevel  = itk.Array[itk.UI]()
shrinkFactorsPerLevel.SetSize(numberOfLevels)
shrinkFactorsPerLevel.Fill(1)

smoothingSigmasPerLevel = itk.Array[itk.UI]()
smoothingSigmasPerLevel.SetSize(numberOfLevels)
smoothingSigmasPerLevel.Fill(0)

registration.SetShrinkFactorsPerLevel(shrinkFactorsPerLevel)
registration.SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel)

print('Initial Value of Metric ', metric.GetValue())
try:
    registration.Update()
except:
    print('Error Occured')

print('Final Value of Metric ', metric.GetValue())

finalTransform = registration.GetModifiableTransform()
numberOfPoints = movingMesh.GetNumberOfPoints()
for n in range(0, numberOfPoints):
    movingMesh.SetPoint(n, finalTransform.TransformPoint(movingMesh.GetPoint(n)))

itk.meshwrite(movingMesh, "synMovingMesh.vtk")
