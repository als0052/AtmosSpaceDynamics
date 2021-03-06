# Orbital Mechanics JSON File

**Last Updated**: 08-31-2020

Details about the ``orbital_mechanics.json`` file

-----

## Purpose

Define variables in JSON syntax to be parsed by the Python program.

## Defining a variable

How the file should be formatted and how a variable should be written out in JSON

```JSON
"variable_name": {
	"render": {
		"latex": "How the variable should be rendered in LaTeX. Must use double back-strokes to properly escape them in JSON.",
		"html": "How the variable should be rendered in HTML",
		"unicode": {
			"oct": "How the variable should be rendered in octal",
			"dec": "How the variable should be rendered in decimal",
			"hex": "How the variable should be rendered in hexadecimal"
		}
	},
	"description": "What does the variable stand for",
	"units": "The units of the variable",
	"definition": "Where in the book(s) is this defined. Provide full citation in APA?",
	"standard_values": "If there are common values associated with this variable they should be listed as name-value pairs"
}
```

### Example definition
```JSON
"mu": {
	"render": {
		"latex": "\\mu",
		"html": "&mu;",
		"unincode": {
			"oct": "01674",
			"dec": "956",
			"hex": "0x3BC"
		}
	},
	"description": "standard gravitational parameter",
	"units": "m^3/s^2",
	"definition": "some page in the book",
	"standard_values": {
		"Sun": 1.32712440018e20,
		"Mercury": 2.2032e13,
		"Venus": 3.24859e14,
		"Earth": 3.986004418e14,
		"moon": 4.9048695e12,
		"Mars": 4.282837e13,
		"Jupiter": 1.26686534e17,
		"Saturn": 3.7931187e16,
		"Uranus": 5.793939e15,
		"Neptune": 6.836529e15,
		"Pluto": 8.71e11,
		"Ceres": 6.26325e10,
		"Eris": 1.108e12
	}
}
```
