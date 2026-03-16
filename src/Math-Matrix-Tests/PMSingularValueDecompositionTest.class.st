"
Please comment me using the following template inspired by Class Responsibility Collaborator (CRC) design:

For the Class part:  State a one line summary. For example, ""I represent a paragraph of text"".

For the Responsibility part: Three sentences about my main responsibilities - what I do, what I know.

For the Collaborators Part: State my main collaborators and one line about how I interact with them. 

Public API and Key Messages
- message one   
- message two 
- (for bonus points) how to create instances.

   One simple example is simply gorgeous.
 
Internal Representation and Key Implementation Points.

    Instance Variables
	ones:		<Object>
	randomMatrix:		<Object>
	s_matrix:		<Object>
	u:		<Object>
	v:		<Object>


    Implementation Points
"
Class {
	#name : 'PMSingularValueDecompositionTest',
	#superclass : 'ParametrizedTestCase',
	#instVars : [
		'example'
	],
	#classInstVars : [
		'example1',
		'example2',
		'example3',
		'example4'
	],
	#category : 'Math-Matrix-Tests',
	#package : 'Math-Matrix-Tests'
}

{ #category : 'initi' }
PMSingularValueDecompositionTest class >> loadExamples [

	example1 := PMSingularValueDecompositionExample
		            newMatrix: (PMMatrix rows: #( #( 0 4 ) #( 0 0 ) ))
		            u: (PMMatrix rows: #( #( 1 0 ) #( 0 1 ) ))
		            v: (PMMatrix rows: #( #( 0 1 ) #( 1 0 ) ))
		            s: (PMMatrix rows: #( #( 4 0 ) #( 0 0 ) )).
	example2 := PMSingularValueDecompositionExample
		            newMatrix: (PMMatrix rows: #( #( 0 1 0 0 ) #( 0 0 2 0 ) #( 0 0 0 3 ) #( 0 0 0 0 ) ))
		            u: (PMMatrix rows: #( #( 0 0 1 0 ) #( 0 1 0 0 ) #( 1 0 0 0 ) #( 0 0 0 1 ) ))
		            v: (PMMatrix rows: #( #( 0 0 0 1 ) #( 0 0 1 0 ) #( 0 1 0 0 ) #( 1 0 0 0 ) ))
		            s: (PMMatrix rows: #( #( 3 0 0 0 ) #( 0 2 0 0 ) #( 0 0 1 0 ) #( 0 0 0 0 ) )).
	example3 := PMSingularValueDecompositionExample
		            newMatrix: (PMMatrix rows: #( #( -1 1 0 ) #( 0 -1 1 ) ))
		            u: (PMMatrix rows: { { (1 / 2 sqrt). (1 / 2 sqrt) }. { (-1 / 2 sqrt). (1 / 2 sqrt) } })
		            v: (PMMatrix rows: { { (-1 / 6 sqrt). (-1 / 2 sqrt). (1 / 3 sqrt) }. { (2 / 6 sqrt). 0. (1 / 3 sqrt) }. { (-1 / 6 sqrt). (1 / 2 sqrt). (1 / 3 sqrt) } })
		            s: (PMMatrix rows: { { 3 sqrt. 0. 0 }. { 0. 1. 0 } }).
	example4 := PMSingularValueDecompositionExample
		            newMatrix: (PMMatrix rows: #( #( 1 0 0 0 2 ) #( 0 0 3 0 0 ) #( 0 0 0 0 0 ) #( 0 2 0 0 0 ) ))
		            u: (PMMatrix rows: #( #( 0 1 0 0 ) #( 1 0 0 0 ) #( 0 0 0 1 ) #( 0 0 1 0 ) ))
		            v: (PMMatrix rows: { { 0. 0.2 sqrt. 0. 0.8 sqrt. 0 }. { 0. 0. 1. 0. 0 }. { 1. 0. 0. 0. 0 }. { 0. 0. 0. 0. 1 }. { 0. 0.8 sqrt. 0. 0.2 sqrt negated. 0 } })
		            s: (PMMatrix rows: { { 3. 0. 0. 0. 0 }. { 0. 5 sqrt. 0. 0. 0 }. { 0. 0. 2. 0. 0 }. { 0. 0. 0. 0. 0 } }).
]

{ #category : 'building suites' }
PMSingularValueDecompositionTest class >> testParameters [

	self loadExamples.

	^ ParametrizedTestMatrix new forSelector: #example addOptions: {
			  example1.
			  example2.
			  example3.
			  example4 }
]

{ #category : 'accessing' }
PMSingularValueDecompositionTest >> actualS [

	^ self example actualS
]

{ #category : 'accessing' }
PMSingularValueDecompositionTest >> actualU [

	^ self example actualU
]

{ #category : 'accessing' }
PMSingularValueDecompositionTest >> actualV [

	^ self example actualV
]

{ #category : 'accessing' }
PMSingularValueDecompositionTest >> example [

	^ example
]

{ #category : 'accessing' }
PMSingularValueDecompositionTest >> example: anObject [

	example := anObject
]

{ #category : 'accessing' }
PMSingularValueDecompositionTest >> matrix [

	^ self example matrix
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testIsExampleCorrect [
	"Verifying that the example SVD is correct"

	| reconstructed |
	reconstructed := self actualU * self actualS * self actualV transpose.
	self assert: reconstructed closeTo: self matrix.
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testIsExampleOrthonormalU [
	"Verifying that matrix U in the given example is orthonormal. That is, all columns of U are unit vectors and orthogonal to each other."
	
	| m identity |
		
	m := self actualU numberOfRows.
	
	identity := PMSymmetricMatrix identity: m.
	self assert: (self actualU transpose * self actualU) closeTo: identity
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testIsExampleOrthonormalV [
	"Verifying that matrix V in the given example is orthonormal. That is, all columns of V are unit vectors and orthogonal to each other."
	
	| n identity |
		
	n := self actualV numberOfRows.
	
	identity := PMSymmetricMatrix identity: n.
	self assert: (self actualV transpose * self actualV) closeTo: identity
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testIsOrthonormalU [

	"Verifying that matrix U is orthonormal. That is, all columns of U are unit vectors and orthogonal to each other."

	| u m identity |
	u := self matrix decomposeSV leftSingularMatrix.
	m := u numberOfRows.
	identity := PMSymmetricMatrix identity: m.
	self assert: u transpose * u closeTo: identity
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testIsOrthonormalV [

	"Verifying that matrix U is orthonormal. That is, all columns of U are unit vectors and orthogonal to each other."

	| v n identity |
	v := self matrix decomposeSV rightSingularMatrix.
	n := v numberOfRows.
	identity := PMSymmetricMatrix identity: n.
	self assert: v transpose * v closeTo: identity
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testIsSquareU [

	"U should be a square matrix"

	| u |
	u := self matrix decomposeSV leftSingularMatrix.
	self assert: u numberOfRows equals: u numberOfColumns
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testIsSquareV [

	"V should be a square matrix"

	| v |
	v := self matrix decomposeSV rightSingularMatrix.
	self assert: v numberOfRows equals: v numberOfColumns
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testReconstruction [

	| svd u v s reconstructed |
	svd := self matrix decomposeSV.
	u := svd leftSingularMatrix.
	v := svd rightSingularMatrix.
	s := svd diagonalSingularValueMatrix.
	reconstructed := u * s * v transpose.
	self assert: reconstructed closeTo: self matrix
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testShapeS [

	"If A is an (m x n) matrix, then S should be an (m x n) matrix"

	| s m n |
	s := self matrix decomposeSV diagonalSingularValueMatrix.
	m := self matrix numberOfRows.
	n := self matrix numberOfColumns.
	self assert: s numberOfRows equals: m.
	self assert: s numberOfColumns equals: n
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testShapeU [

	"If A is an (m x n) matrix, then U should be an (m x m) matrix"

	| u m |
	u := self matrix decomposeSV leftSingularMatrix.
	m := self matrix numberOfRows.
	self assert: u numberOfRows equals: m.
	self assert: u numberOfColumns equals: m
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testShapeV [

	"If A is an (m x n) matrix, then V should be an (n x n) matrix"

	| v n |
	v := self matrix decomposeSV rightSingularMatrix.
	n := self matrix numberOfColumns.
	self assert: v numberOfRows equals: n.
	self assert: v numberOfColumns equals: n
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testValueS [

	"Comparing S to its known value"

	| s |
	s := self matrix decomposeSV diagonalSingularValueMatrix.
	self assert: s closeTo: self actualS
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testValueU [

	"Comparing U to its known value"

	| u |
	u := self matrix decomposeSV leftSingularMatrix.
	self assert: u closeTo: self actualU
]

{ #category : 'tests' }
PMSingularValueDecompositionTest >> testValueV [

	"Comparing V to its known value"

	| v |
	v := self matrix decomposeSV rightSingularMatrix.
	self assert: v closeTo: self actualV
]
