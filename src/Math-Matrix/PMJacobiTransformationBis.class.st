Class {
	#name : 'PMJacobiTransformationBis',
	#superclass : 'PMJacobiTransformation',
	#instVars : [
		'n',
		'd',
		'v',
		'nrot',
		'nMax'
	],
	#category : 'Math-Matrix',
	#package : 'Math-Matrix'
}

{ #category : 'as yet unclassified' }
PMJacobiTransformationBis >> eigsrt [

	| k p |
	1 to: n - 1 do: [ :i |
			k := i.
			p := d at: i.
			i + 1 to: n do: [ :j |
					(d at: j) >= p ifTrue: [
							k := j.
							p := d at: j ] ].
			k ~= i ifTrue: [
					d at: k put: (d at: i).
					d at: i put: p.
					1 to: n do: [ :j |
							p := v at: j at: i.
							v at: j at: i put: (v at: j at: k).
							v at: j at: k put: p ] ] ]
]

{ #category : 'operation' }
PMJacobiTransformationBis >> evaluate [

	^ d
]

{ #category : 'initialization' }
PMJacobiTransformationBis >> initialize: aMatrix [

	| a b z sm tresh c s g h t theta tau |
	n := aMatrix numberOfColumns.
	a := PMMatrix new: n.
	v := PMMatrix zerosRows: n cols: n.
	1 to: n do: [ :i |
			1 to: n do: [ :j |
					a at: i at: j put: (aMatrix at: i at: j).
					i = j ifTrue: [ v at: i at: j put: 1 ] ] ].
	nrot := 0.
	nMax := 500.

	b := PMVector new: n.
	d := PMVector new: n.

	1 to: n do: [ :ip |
			b at: ip put: (a at: ip at: ip).
			d at: ip put: (b at: ip) ].
	z := PMVector zeros: n.

	1 to: 50 do: [ :i |
			sm := self sumUpperTriangle: a.
			sm = 0 ifTrue: [
					self eigsrt.
					^ self ].
			tresh := i < 4
				         ifTrue: [ 0.2 * sm / n squared ]
				         ifFalse: [ 0 ].
			1 to: n - 1 do: [ :ip |
					ip + 1 to: n do: [ :iq |
							g := 100 * (a at: ip at: iq) abs.
							i > 4 & ((d at: ip) abs + g = (d at: ip) abs) & ((d at: iq) abs + g = (d at: iq) abs)
								ifTrue: [ a at: ip at: iq put: 0 ]
								ifFalse: [
										(a at: ip at: iq) abs > tresh ifTrue: [
												h := (d at: iq) - (d at: ip).
												h abs + g = h abs
													ifTrue: [ t := (a at: ip at: iq) / h ]
													ifFalse: [
															theta := 0.5 * h / (a at: ip at: iq).
															t := 1 / (theta abs + (1.0 + theta squared) sqrt).
															theta < 0 ifTrue: [ t := t negated ] ].
												c := 1 / (1 + t squared) sqrt.
												s := t * c.
												tau := s / (1 + c).
												h := t * (a at: ip at: iq).
												z at: ip put: (z at: ip) - h.
												z at: iq put: (z at: iq) + h.
												d at: ip put: (d at: ip) - h.
												d at: iq put: (d at: iq) + h.
												a at: ip at: iq put: 0.
												1 to: ip - 1 do: [ :j |
														g := a at: j at: ip.
														h := a at: j at: iq.
														a at: j at: ip put: g - (s * (h + (g * tau))).
														a at: j at: iq put: h + (s * (g - (h * tau))) ].
												ip + 1 to: iq - 1 do: [ :j |
														g := a at: ip at: j.
														h := a at: j at: iq.
														a at: ip at: j put: g - (s * (h + (g * tau))).
														a at: j at: iq put: h + (s * (g - (h * tau))) ].
												iq + 1 to: n do: [ :j |
														g := a at: ip at: j.
														h := a at: iq at: j.
														a at: ip at: j put: g - (s * (h + (g * tau))).
														a at: iq at: j put: h + (s * (g - (h * tau))) ].
												1 to: n do: [ :j |
														g := v at: j at: ip.
														h := v at: j at: iq.
														v at: j at: ip put: g - (s * (h + (g * tau))).
														v at: j at: iq put: h + (s * (g - (h * tau))) ].
												nrot := nrot + 1 ] ] ] ].
			1 to: n do: [ :ip |
					b at: ip put: (b at: ip) + (z at: ip).
					d at: ip put: (b at: ip).
					z at: ip put: 0 ] ].
	self eigsrt.
	^ self
]

{ #category : 'as yet unclassified' }
PMJacobiTransformationBis >> sumUpperTriangle: a [

	| sm |
	sm := 0.
	1 to: n - 1 do: [ :ip | ip + 1 to: n do: [ :iq | sm := sm + (a at: ip at: iq) abs ] ].
	^ sm
]

{ #category : 'private - transforming' }
PMJacobiTransformationBis >> transform [

	^ v
]
