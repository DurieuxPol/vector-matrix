Class {
	#name : 'PMSingularValueDecompositionExample',
	#superclass : 'Object',
	#instVars : [
		'matrix',
		'actualU',
		'actualV',
		'actualS'
	],
	#category : 'Math-Matrix-Tests',
	#package : 'Math-Matrix-Tests'
}

{ #category : 'instance creation' }
PMSingularValueDecompositionExample class >> newMatrix: aMatrix u: uMatrix v: vMatrix s: sMatrix [

	^ self new
		  matrix: aMatrix;
		  actualU: uMatrix;
		  actualV: vMatrix;
		  actualS: sMatrix;
		  yourself
]

{ #category : 'accessing' }
PMSingularValueDecompositionExample >> actualS [

	^ actualS
]

{ #category : 'accessing' }
PMSingularValueDecompositionExample >> actualS: anObject [

	actualS := anObject
]

{ #category : 'accessing' }
PMSingularValueDecompositionExample >> actualU [

	^ actualU
]

{ #category : 'accessing' }
PMSingularValueDecompositionExample >> actualU: anObject [

	actualU := anObject
]

{ #category : 'accessing' }
PMSingularValueDecompositionExample >> actualV [

	^ actualV
]

{ #category : 'accessing' }
PMSingularValueDecompositionExample >> actualV: anObject [

	actualV := anObject
]

{ #category : 'accessing' }
PMSingularValueDecompositionExample >> matrix [

	^ matrix
]

{ #category : 'accessing' }
PMSingularValueDecompositionExample >> matrix: anObject [

	matrix := anObject
]
