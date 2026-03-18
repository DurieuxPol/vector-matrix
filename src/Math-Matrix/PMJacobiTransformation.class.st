"
Based on the algorithm in ""Numerical Recipes in Fortran 77: The Art of Scientific Computing"", Chapter 11, page 460 ( https://perso.ens-lyon.fr/christophe.winisdoerffer/INTRO_NUM/NumericalRecipesinF77.pdf ).


"
Class {
	#name : 'PMJacobiTransformation',
	#superclass : 'Object',
	#instVars : [
		'numberOfRotations',
		'matrixSize',
		'eigenValues',
		'eigenVectors'
	],
	#category : 'Math-Matrix',
	#package : 'Math-Matrix'
}

{ #category : 'creation' }
PMJacobiTransformation class >> matrix: aSymmetricMatrix [

	^ super new initialize: aSymmetricMatrix
]

{ #category : 'accessing' }
PMJacobiTransformation >> eigenValues [

	^ eigenValues
]

{ #category : 'accessing' }
PMJacobiTransformation >> eigenVectors [

	^ eigenVectors
]

{ #category : 'computing' }
PMJacobiTransformation >> evaluateIterationsOn: jacobiMatrix accumulations: accumulationVector corrections: correctionVector [

	| sum |
	1 to: 50 do: [ :i |
			sum := self sumUpperTriangle: jacobiMatrix.
			sum = 0 ifTrue: [ ^ self ].
			self
				performJacobiIterationOn: jacobiMatrix
				corrections: correctionVector
				sum: sum
				iteration: i.
			1 to: matrixSize do: [ :ip |
					accumulationVector at: ip put: (accumulationVector at: ip) + (correctionVector at: ip).
					eigenValues at: ip put: (accumulationVector at: ip).
					correctionVector at: ip put: 0 ] ]
]

{ #category : 'initialization' }
PMJacobiTransformation >> initialize: inputMatrix [

	| jacobiMatrix accumulationVector correctionVector |
	matrixSize := inputMatrix numberOfColumns.
	jacobiMatrix := PMMatrix new: matrixSize.
	eigenVectors := PMMatrix zerosRows: matrixSize cols: matrixSize.
	1 to: matrixSize do: [ :i |
			1 to: matrixSize do: [ :j |
					jacobiMatrix at: i at: j put: (inputMatrix at: i at: j).
					i = j ifTrue: [ eigenVectors at: i at: j put: 1 ] ] ].
	numberOfRotations := 0.

	accumulationVector := PMVector new: matrixSize.
	eigenValues := Array new: matrixSize.

	1 to: matrixSize do: [ :ip |
			accumulationVector at: ip put: (jacobiMatrix at: ip at: ip).
			eigenValues at: ip put: (accumulationVector at: ip) ].

	correctionVector := PMVector zeros: matrixSize.

	self evaluateIterationsOn: jacobiMatrix accumulations: accumulationVector corrections: correctionVector.
	self sortEigen.

	^ self
]

{ #category : 'computing' }
PMJacobiTransformation >> performJacobiIterationOn: jacobiMatrix corrections: correctionVector sum: sum iteration: iter [

	| tresh |
	tresh := iter < 4
		         ifTrue: [ 0.2 * sum / matrixSize squared ]
		         ifFalse: [ 0 ].
	1 to: matrixSize - 1 do: [ :ip |
			ip + 1 to: matrixSize do: [ :iq |
					| c s g h t theta tau |
					g := 100 * (jacobiMatrix at: ip at: iq) abs.
					iter > 4 & ((eigenValues at: ip) abs + g = (eigenValues at: ip) abs)
					& ((eigenValues at: iq) abs + g = (eigenValues at: iq) abs)
						ifTrue: [ jacobiMatrix at: ip at: iq put: 0 ]
						ifFalse: [
								(jacobiMatrix at: ip at: iq) abs > tresh ifTrue: [
										h := (eigenValues at: iq) - (eigenValues at: ip).
										h abs + g = h abs
											ifTrue: [ t := (jacobiMatrix at: ip at: iq) / h ]
											ifFalse: [
													theta := 0.5 * h / (jacobiMatrix at: ip at: iq).
													t := 1 / (theta abs + (1.0 + theta squared) sqrt).
													theta < 0 ifTrue: [ t := t negated ] ].
										c := 1 / (1 + t squared) sqrt.
										s := t * c.
										tau := s / (1 + c).
										h := t * (jacobiMatrix at: ip at: iq).
										correctionVector at: ip put: (correctionVector at: ip) - h.
										correctionVector at: iq put: (correctionVector at: iq) + h.
										eigenValues at: ip put: (eigenValues at: ip) - h.
										eigenValues at: iq put: (eigenValues at: iq) + h.
										jacobiMatrix at: ip at: iq put: 0.
										self
											performRotationIn: jacobiMatrix
											atIndexes: ip
											and: iq
											using: s
											and: tau ] ] ] ]
]

{ #category : 'computing' }
PMJacobiTransformation >> performRotationIn: jacobiMatrix atIndexes: ip and: iq using: s and: tau [

	| g h |
	1 to: ip - 1 do: [ :j |
			g := jacobiMatrix at: j at: ip.
			h := jacobiMatrix at: j at: iq.
			jacobiMatrix at: j at: ip put: g - (s * (h + (g * tau))).
			jacobiMatrix at: j at: iq put: h + (s * (g - (h * tau))) ].
	ip + 1 to: iq - 1 do: [ :j |
			g := jacobiMatrix at: ip at: j.
			h := jacobiMatrix at: j at: iq.
			jacobiMatrix at: ip at: j put: g - (s * (h + (g * tau))).
			jacobiMatrix at: j at: iq put: h + (s * (g - (h * tau))) ].
	iq + 1 to: matrixSize do: [ :j |
			g := jacobiMatrix at: ip at: j.
			h := jacobiMatrix at: iq at: j.
			jacobiMatrix at: ip at: j put: g - (s * (h + (g * tau))).
			jacobiMatrix at: iq at: j put: h + (s * (g - (h * tau))) ].
	1 to: matrixSize do: [ :j |
			g := eigenVectors at: j at: ip.
			h := eigenVectors at: j at: iq.
			eigenVectors at: j at: ip put: g - (s * (h + (g * tau))).
			eigenVectors at: j at: iq put: h + (s * (g - (h * tau))) ].
	numberOfRotations := numberOfRotations + 1
]

{ #category : 'sorting' }
PMJacobiTransformation >> sortEigen [
	"Rounds the eigen values close to 0, and sorts the eigen values (and their respective eigen vectors) in descending order"

	| k value |
	eigenValues := eigenValues collect: [ :each |
			               (each closeTo: 0)
				               ifTrue: [ 0 ]
				               ifFalse: [ each ] ].
	1 to: matrixSize - 1 do: [ :i |
			k := i.
			value := eigenValues at: i.
			i + 1 to: matrixSize do: [ :j |
					(eigenValues at: j) > value ifTrue: [
							k := j.
							value := eigenValues at: j ] ].
			k ~= i ifTrue: [
					eigenValues at: k put: (eigenValues at: i).
					eigenValues at: i put: value.
					1 to: matrixSize do: [ :j |
							value := eigenVectors at: j at: i.
							eigenVectors at: j at: i put: (eigenVectors at: j at: k).
							eigenVectors at: j at: k put: value ] ] ]
]

{ #category : 'computing' }
PMJacobiTransformation >> sumUpperTriangle: a [

	| sm |
	sm := 0.
	1 to: matrixSize - 1 do: [ :ip | ip + 1 to: matrixSize do: [ :iq | sm := sm + (a at: ip at: iq) abs ] ].
	^ sm
]
