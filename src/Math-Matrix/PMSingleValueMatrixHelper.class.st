Class {
	#name : 'PMSingleValueMatrixHelper',
	#superclass : 'PMJacobiTransformationHelper',
	#category : 'Math-Matrix',
	#package : 'Math-Matrix'
}

{ #category : 'initialization' }
PMSingleValueMatrixHelper >> initialize: aSingleValueMatrix [

	eigenvalues := Array with: (aSingleValueMatrix at: 1 at: 1).
	eigenvectors := PMVector with: (PMVector ones: 1)
]
