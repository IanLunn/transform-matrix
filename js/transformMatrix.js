/**
 *
 */
var transformMatrix = (function(Modernizr, element) {

  var self = {};

  /**
   * Convert a CSS transformMatrix function to an array.
   * "matrix(-1, 0, 0, -1, 100, 100)" -> ["-1", "0", "0", "-1", "100", "100"]
   *
   * @param {String} cssTransformMatrix - The raw CSS matrix function
   * @return {Array} transformMatrix - An array representing the matrix
   */
  self.convertTransformMatrix = function(cssTransformMatrix) {

    var transformMatrix = 0;

    // Get the matrix array if the CSS property is a matrix and not "none"
    // for example
    if(cssTransformMatrix.indexOf("matrix") > -1) {

      // Strip the "matrix()" portion of the CSS transform matrix, leaving
      // just the values, then split them into an array
      transformMatrix = cssTransformMatrix.split('(')[1];
      transformMatrix = transformMatrix.split(')')[0];
      transformMatrix = transformMatrix.split(', ');

      var transformMatrixLength = transformMatrix.length;

      // Convert the matrix to numbers
      for(var i = 0; i < transformMatrixLength; i++) {
        var value = transformMatrix[i];

        transformMatrix[i] = parseFloat(value);
      }
    }

    return transformMatrix;
  };

  /**
   *
   */
  self.determinant = function(matrix) {

    return matrix.m14 * matrix.m23 * matrix.m32 * matrix.m41-matrix.m13 * matrix.m24 * matrix.m32 * matrix.m41 -
matrix.m14 * matrix.m22 * matrix.m33 * matrix.m41+matrix.m12 * matrix.m24 * matrix.m33 * matrix.m41 +
matrix.m13 * matrix.m22 * matrix.m34 * matrix.m41-matrix.m12 * matrix.m23 * matrix.m34 * matrix.m41 -
matrix.m14 * matrix.m23 * matrix.m31 * matrix.m42+matrix.m13 * matrix.m24 * matrix.m31 * matrix.m42 +
matrix.m14 * matrix.m21 * matrix.m33 * matrix.m42-matrix.m11 * matrix.m24 * matrix.m33 * matrix.m42 -
matrix.m13 * matrix.m21 * matrix.m34 * matrix.m42+matrix.m11 * matrix.m23 * matrix.m34 * matrix.m42 +
matrix.m14 * matrix.m22 * matrix.m31 * matrix.m43-matrix.m12 * matrix.m24 * matrix.m31 * matrix.m43 -
matrix.m14 * matrix.m21 * matrix.m32 * matrix.m43+matrix.m11 * matrix.m24 * matrix.m32 * matrix.m43 +
matrix.m12 * matrix.m21 * matrix.m34 * matrix.m43-matrix.m11 * matrix.m22 * matrix.m34 * matrix.m43 -
matrix.m13 * matrix.m22 * matrix.m31 * matrix.m44+matrix.m12 * matrix.m23 * matrix.m31 * matrix.m44 +
matrix.m13 * matrix.m21 * matrix.m32 * matrix.m44-matrix.m11 * matrix.m23 * matrix.m32 * matrix.m44 -
matrix.m12 * matrix.m21 * matrix.m33 * matrix.m44+matrix.m11 * matrix.m22 * matrix.m33 * matrix.m44;
  };

  /**
   *
   */
  self.transpose = function(matrix) {

    for (var n = 0; n <= 4-2; n++) {
      for (var m = n + 1; m <= 4-1; m++) {
        matrix['m'+(n+1)+(m+1)] = matrix['m'+(m+1)+(n+1)];
        matrix['m'+(m+1)+(n+1)] = matrix['m'+(n+1)+(m+1)];
      }
    }

    return matrix;
  };

  /**
   *
   */
  self.transformVector = function(matrix, vector) {

    return matrix.m11* vector.x + matrix.m12* vector.y + matrix.m13* vector.z,
           matrix.m21* vector.x + matrix.m22* vector.y + matrix.m23* vector.z,
           matrix.m31* vector.x + matrix.m32* vector.y + matrix.m33* vector.z;
  };

  /**
   * Ensure that values are not undefined
   */

  self.checkValues = function(row) {
    row.x = row.x ? row.x : 0;
    row.y = row.y ? row.y : 0;
    row.z = row.z ? row.z : 0;
    row.w = row.w ? row.w : 0;
  };

  /**
   *
   */
  self.length = function(row) {

    	this.checkValues(row);

      return Math.sqrt(row.x * row.x + row.y * row.y + row.z * row.z);
  };

  /**
   * Get a normalised representation of the vector
   */

  self.normalise = function(row) {
    var len = this.length(row);

    return {
      "x": row.x / len,
      "y": row.y / len,
      "z": row.z / len,
      "w": 0
    };
  };

  /**
   *
   */
  self.dot = function(row0, row1) {

    return row0.x * row1.x + row0.y * row1.y + row0.z * row1.z + row0.w * row1.w;
  };

  self.cross = function(row1, row2) {

		return {
      "x": row1.y*row2.z - row1.z*row2.y,
      "y": row1.z*row2.x - row1.x*row2.z,
      "z": row1.x*row2.y - row1.y*row2.x,
      "w": 0
    };
	};

  /**
   *
   */
  self.combine = function(row1, aPoint, ascl, bscl) {

    return {
      "x": (ascl * row1.x) + (bscl * aPoint.x),
      "y": (ascl * row1.y) + (bscl * aPoint.y),
      "z": (ascl * row1.z) + (bscl * aPoint.z),
      "w": 0
    };
  };

  self.goto = function(id) {

    // zero-base
    var stepId = id - 1,
        currentStep = self.steps[stepId];

    var styles = getComputedStyle(currentStep, null)

    // Get the transform matrix of the step
    // Convert the transformMatrix to an array
    var cssTransformMatrix = styles[Modernizr.prefixed("transform")],
        transformMatrix = self.convertTransformMatrix(cssTransformMatrix),
        transformMatrixLength = transformMatrix.length;

    var translateX = 0,
        translateY = 0,
        translateZ = 0,
        scaleX = 1,
        scaleY = 1,
        scaleZ = 1,
        skewX = 0,
        skewY = 0,
        skewZ = 0,
        angle = 0,
        rotateX = 0,
        rotateY = 0,
        rotateZ = 0,
        stepX = 0,
        stepY = 0,
        stepZ = 0,
        quaternion = {};

    // Get the translate X and Y positions of the next step and them to
    // the offsetLeft and offsetTop properties
    if(transformMatrix !== 0) {

      // 2D Matrix
      if (transformMatrix.length === 6) {

        var row0x = transformMatrix[0];
        var row0y = transformMatrix[1];
        var row1x = transformMatrix[2];
        var row1y = transformMatrix[3];

        translateX = transformMatrix[4];
        translateY = transformMatrix[5];

        scaleX = Math.sqrt(row0x * row0x + row0y * row0y);
        scaleY = Math.sqrt(row1x * row1x + row1y * row1y);

        // If determinant is negative, one axis was flipped
        var determinant = row0x * row1y - row0y * row1x;
        if (determinant < 0) {

          // Flip axis with minimum unit vector dot product
          if (row0x < row1y) {
            scaleX = -scaleX;
          } else {
            scaleY = -scaleY;
          }
        }

        // Compute rotation and renormalize matrix
        angle = Math.atan2(row0y, row0x);
      }

      // 3D Matrix
      else if (transformMatrix.length === 16) {

        var matrix = {};
        matrix.m11 = transformMatrix[0];
        matrix.m12 = transformMatrix[1];
        matrix.m13 = transformMatrix[2];
        matrix.m14 = transformMatrix[3];
        matrix.m21 = transformMatrix[4];
        matrix.m22 = transformMatrix[5];
        matrix.m23 = transformMatrix[6];
        matrix.m24 = transformMatrix[7];
        matrix.m31 = transformMatrix[8];
        matrix.m32 = transformMatrix[9];
        matrix.m33 = transformMatrix[10];
        matrix.m34 = transformMatrix[11];
        matrix.m41 = transformMatrix[12];
        matrix.m42 = transformMatrix[13];
        matrix.m43 = transformMatrix[14];
        matrix.m44 = transformMatrix[15];

        // Normalize the matrix.
        if (matrix.m44 === 0) {
          return false;
        }

        for (var i = 1; i <= 4; i++) {
          for (var j = 1; j <= 4; j++) {
            matrix['m'+i+j] /= matrix.m44;
          }
        }

        // perspectiveMatrix is used to solve for perspective, but it also
        // provides an easy way to test for singularity of the upper 3x3
        // component.
        var perspectiveMatrix = transformMatrix;

        for (i = 1; i <= 3; i++) {
          perspectiveMatrix['m'+i+'4'] = 0;
        }

        perspectiveMatrix.m44 = 1;

        if (self.determinant(perspectiveMatrix) === 0) {
          return false;
        }

        // First, isolate perspective
        if (matrix.m14 !== 0 || matrix.m24 !== 0 || matrix.m34 !== 0) {

          var rightHandSide = [];

          // rightHandSide is the right hand side of the equation
          rightHandSide[0] = matrix.m14;
          rightHandSide[1] = matrix.m24;
          rightHandSide[2] = matrix.m34;
          rightHandSide[3] = matrix.m44;

          // Solve the equation by inverting perspectiveMatrix and
          // multiplying rightHandSide by the inverse
          var inversePerspectiveMatrix = self.inverse(perspectiveMatrix);
          var transposedInversePerspectiveMatrix = self.transpose(inversePerspectiveMatrix);
          var perspective = transformVector(transposedInversePerspectiveMatrix, rightHandSide);

          // Clear the perspective partition
          matrix.m14 = 0;
          matrix.m24 = 0;
          matrix.m34 = 0;
          matrix.m44 = 1;
        }

        //  No perspective
        else {

          var perspective = [];

          perspective[0] = 0;
          perspective[1] = 0;
          perspective[2] = 0;
          perspective[3] = 1;
        }

        // Next take care of translation
        translateX = matrix.m41;
        translateY = matrix.m42;
        translateZ = matrix.m43;

        matrix.m41 = 0;
        matrix.m42 = 0;
        matrix.m43 = 0;

        var row = [];

        // Now get scale and shear. 'row' is a 3 element array of 3 component vectors
        for (i = 1; i <= 3; i++) {

          row[i - 1] = {};
          row[i - 1].x = matrix['m'+i+'1'];
          row[i - 1].y = matrix['m'+i+'2'];
          row[i - 1].z = matrix['m'+i+'3'];
          row[i - 1].w = 0;
        }

        // Compute X scale factor and normalize first row.
        scaleX = self.length(row[0]);
        row[0] = self.normalise(row[0]);

        // Compute XY shear factor and make 2nd row orthogonal to 1st.
        skewX = self.dot(row[0], row[1]);
        row[1] = self.combine(row[1], row[0], 1.0, -skewX);

        // Now, compute Y scale and normalize 2nd row.
        scaleY = self.length(row[1]);
        row[1] = self.normalise(row[1]);
        skewX /= scaleY;

        // Compute XZ and YZ shears, orthogonalize 3rd row
        skewY = self.dot(row[0], row[2]);

        row[2] = self.combine(row[2], row[0], 1.0, -skewY);
        skewZ = self.dot(row[1], row[2]);

        row[2] = self.combine(row[2], row[1], 1.0, -skewZ);

        // Next, get Z scale and normalize 3rd row.
        scaleZ = self.length(row[2]);
        row[2] = self.normalise(row[2]);
        skewY /= scaleZ;

        // At this point, the matrix (in rows) is orthonormal.
        // Check for a coordinate system flip.  If the determinant
        // is -1, then negate the matrix and the scaling factors.
        var pdum3 = self.cross(row[1], row[2]);

        if (self.dot(row[0], pdum3) < 0) {
          for (i = 0; i < 3; i++) {
            if(i === 0) {
              scaleX *= -1;
            }

            if(i === 1) {
              scaleY *= -1;
            }

            if(i === 2) {
              scaleZ *= -1;
            }

            row[i].x *= -1;
            row[i].y *= -1;
            row[i].z *= -1;
          }
        }

        // Now, get the rotations out
        quaternion.x = 0.5 * Math.sqrt(Math.max(1 + row[0].x - row[1].y - row[2].z, 0));
        quaternion.y = 0.5 * Math.sqrt(Math.max(1 - row[0].x + row[1].y - row[2].z, 0));
        quaternion.z = 0.5 * Math.sqrt(Math.max(1 - row[0].x - row[1].y + row[2].z, 0));
        quaternion.w = 0.5 * Math.sqrt(Math.max(1 + row[0].x + row[1].y + row[2].z, 0));

        if (row[2].y > row[1].z) {
          quaternion.x = -quaternion.x;
        }

        if (row[0].z > row[2].x) {
          quaternion.y = -quaternion.y;
        }

        if (row[1].x > row[0].y) {
          quaternion.z = -quaternion.z;
        }

        var scaleMatrix = $M([
          [1,0,0,0],
          [0,1,0,0],
          [0,0,1,0],
          [0,0,0,1]
        ]);

        var rotateXMatrix = $M([
          [1,0,0,0],
          [0,Math.cos(quaternion.x), Math.sin(-quaternion.x), 0],
          [0,Math.sin(quaternion.x), Math.cos( quaternion.x), 0],
          [0,0,0,1]
        ]);

        var rotateYMatrix = $M([
          [Math.cos(quaternion.y), 0, Math.sin(quaternion.y),0],
          [0,1,0,0],
          [Math.sin(-quaternion.y), 0, Math.cos(quaternion.y), 0],
          [0,0,0,1]
        ]);

        var rotateZMatrix = $M([
          [Math.cos(quaternion.z), Math.sin(-quaternion.z), 0, 0],
          [Math.sin(quaternion.z), Math.cos( quaternion.z), 0, 0],
          [0,0,1,0],
          [0,0,0,1]
        ]);

        var translateMatrix = $M([
          [1,0,0,0],
          [0,1,0,0],
          [0,0,1,0],
          [-translateX,-translateY,-translateZ,1]
        ]);

        var transformMatrix = rotateXMatrix.x(rotateYMatrix).x(rotateZMatrix).x(scaleMatrix).x(translateMatrix);

        var transformStyle = "matrix3d(";
            transformStyle += transformMatrix.e(1,1).toFixed(10) + ","
                        + transformMatrix.e(1,2).toFixed(10) + ","
                        + transformMatrix.e(1,3).toFixed(10) + ","
                        + transformMatrix.e(1,4).toFixed(10) + ",";

            transformStyle += transformMatrix.e(2,1).toFixed(10) + ","
                        + transformMatrix.e(2,2).toFixed(10) + ","
                        + transformMatrix.e(2,3).toFixed(10) + ","
                        + transformMatrix.e(2,4).toFixed(10) + ",";

            transformStyle += transformMatrix.e(3,1).toFixed(10) + ","
                        + transformMatrix.e(3,2).toFixed(10) + ","
                        + transformMatrix.e(3,3).toFixed(10) + ","
                        + transformMatrix.e(3,4).toFixed(10) + ",";

            transformStyle += transformMatrix.e(4,1).toFixed(10) + ","
                        + transformMatrix.e(4,2).toFixed(10) + ","
                        + transformMatrix.e(4,3).toFixed(10) + ","
                        + transformMatrix.e(4,4).toFixed(10);

            transformStyle += ")";
      }
    }

    canvas.style[Modernizr.prefixed("transitionDuration")] = "500ms";
    canvas.style[Modernizr.prefixed("transformOrigin")] = "left top 0";
    canvas.style[Modernizr.prefixed("transformProperty")] = Modernizr.prefixed("transform") + ", " + Modernizr.prefixed("transformOrigin");
    canvas.style[Modernizr.prefixed("transform")] = transformStyle;

    // Cycle between steps
    if(self.currentStepId !== self.noOfSteps) {
      self.currentStepId++;
    }else{
      self.currentStepId = 1;
    }
  };

  var canvas = document.getElementById("canvas"),
      nextButton = document.getElementById("next");

  self.steps = canvas.querySelectorAll(".step");
  self.noOfSteps = self.steps.length;
  self.currentStepId = 1;


  nextButton.addEventListener("click", function() {
    self.goto(self.currentStepId);
  });

  self.goto(self.currentStepId);
});
