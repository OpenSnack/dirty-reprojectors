// node_modules/d3-array/src/fsum.js
var Adder = class {
  constructor() {
    this._partials = new Float64Array(32);
    this._n = 0;
  }
  add(x) {
    const p = this._partials;
    let i = 0;
    for (let j = 0; j < this._n && j < 32; j++) {
      const y = p[j], hi = x + y, lo = Math.abs(x) < Math.abs(y) ? x - (hi - y) : y - (hi - x);
      if (lo)
        p[i++] = lo;
      x = hi;
    }
    p[i] = x;
    this._n = i + 1;
    return this;
  }
  valueOf() {
    const p = this._partials;
    let n = this._n, x, y, lo, hi = 0;
    if (n > 0) {
      hi = p[--n];
      while (n > 0) {
        x = hi;
        y = p[--n];
        hi = x + y;
        lo = y - (hi - x);
        if (lo)
          break;
      }
      if (n > 0 && (lo < 0 && p[n - 1] < 0 || lo > 0 && p[n - 1] > 0)) {
        y = lo * 2;
        x = hi + y;
        if (y == x - hi)
          hi = x;
      }
    }
    return hi;
  }
};

// node_modules/d3-array/src/merge.js
function* flatten(arrays) {
  for (const array of arrays) {
    yield* array;
  }
}
function merge(arrays) {
  return Array.from(flatten(arrays));
}

// node_modules/d3-array/src/range.js
function range(start, stop, step) {
  start = +start, stop = +stop, step = (n = arguments.length) < 2 ? (stop = start, start = 0, 1) : n < 3 ? 1 : +step;
  var i = -1, n = Math.max(0, Math.ceil((stop - start) / step)) | 0, range3 = new Array(n);
  while (++i < n) {
    range3[i] = start + i * step;
  }
  return range3;
}

// node_modules/d3-geo/src/math.js
var epsilon = 1e-6;
var epsilon2 = 1e-12;
var pi = Math.PI;
var halfPi = pi / 2;
var quarterPi = pi / 4;
var tau = pi * 2;
var degrees = 180 / pi;
var radians = pi / 180;
var abs = Math.abs;
var atan = Math.atan;
var atan2 = Math.atan2;
var cos = Math.cos;
var exp = Math.exp;
var hypot = Math.hypot;
var log = Math.log;
var pow = Math.pow;
var sin = Math.sin;
var sign = Math.sign || function(x) {
  return x > 0 ? 1 : x < 0 ? -1 : 0;
};
var sqrt = Math.sqrt;
var tan = Math.tan;
function acos(x) {
  return x > 1 ? 0 : x < -1 ? pi : Math.acos(x);
}
function asin(x) {
  return x > 1 ? halfPi : x < -1 ? -halfPi : Math.asin(x);
}
function haversin(x) {
  return (x = sin(x / 2)) * x;
}

// node_modules/d3-geo/src/noop.js
function noop() {
}

// node_modules/d3-geo/src/stream.js
function streamGeometry(geometry, stream) {
  if (geometry && streamGeometryType.hasOwnProperty(geometry.type)) {
    streamGeometryType[geometry.type](geometry, stream);
  }
}
var streamObjectType = {
  Feature: function(object, stream) {
    streamGeometry(object.geometry, stream);
  },
  FeatureCollection: function(object, stream) {
    var features = object.features, i = -1, n = features.length;
    while (++i < n)
      streamGeometry(features[i].geometry, stream);
  }
};
var streamGeometryType = {
  Sphere: function(object, stream) {
    stream.sphere();
  },
  Point: function(object, stream) {
    object = object.coordinates;
    stream.point(object[0], object[1], object[2]);
  },
  MultiPoint: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n)
      object = coordinates[i], stream.point(object[0], object[1], object[2]);
  },
  LineString: function(object, stream) {
    streamLine(object.coordinates, stream, 0);
  },
  MultiLineString: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n)
      streamLine(coordinates[i], stream, 0);
  },
  Polygon: function(object, stream) {
    streamPolygon(object.coordinates, stream);
  },
  MultiPolygon: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n)
      streamPolygon(coordinates[i], stream);
  },
  GeometryCollection: function(object, stream) {
    var geometries = object.geometries, i = -1, n = geometries.length;
    while (++i < n)
      streamGeometry(geometries[i], stream);
  }
};
function streamLine(coordinates, stream, closed) {
  var i = -1, n = coordinates.length - closed, coordinate;
  stream.lineStart();
  while (++i < n)
    coordinate = coordinates[i], stream.point(coordinate[0], coordinate[1], coordinate[2]);
  stream.lineEnd();
}
function streamPolygon(coordinates, stream) {
  var i = -1, n = coordinates.length;
  stream.polygonStart();
  while (++i < n)
    streamLine(coordinates[i], stream, 1);
  stream.polygonEnd();
}
function stream_default(object, stream) {
  if (object && streamObjectType.hasOwnProperty(object.type)) {
    streamObjectType[object.type](object, stream);
  } else {
    streamGeometry(object, stream);
  }
}

// node_modules/d3-geo/src/area.js
var areaRingSum = new Adder();
var areaSum = new Adder();
var lambda00;
var phi00;
var lambda0;
var cosPhi0;
var sinPhi0;
var areaStream = {
  point: noop,
  lineStart: noop,
  lineEnd: noop,
  polygonStart: function() {
    areaRingSum = new Adder();
    areaStream.lineStart = areaRingStart;
    areaStream.lineEnd = areaRingEnd;
  },
  polygonEnd: function() {
    var areaRing = +areaRingSum;
    areaSum.add(areaRing < 0 ? tau + areaRing : areaRing);
    this.lineStart = this.lineEnd = this.point = noop;
  },
  sphere: function() {
    areaSum.add(tau);
  }
};
function areaRingStart() {
  areaStream.point = areaPointFirst;
}
function areaRingEnd() {
  areaPoint(lambda00, phi00);
}
function areaPointFirst(lambda, phi) {
  areaStream.point = areaPoint;
  lambda00 = lambda, phi00 = phi;
  lambda *= radians, phi *= radians;
  lambda0 = lambda, cosPhi0 = cos(phi = phi / 2 + quarterPi), sinPhi0 = sin(phi);
}
function areaPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  phi = phi / 2 + quarterPi;
  var dLambda = lambda - lambda0, sdLambda = dLambda >= 0 ? 1 : -1, adLambda = sdLambda * dLambda, cosPhi = cos(phi), sinPhi = sin(phi), k2 = sinPhi0 * sinPhi, u = cosPhi0 * cosPhi + k2 * cos(adLambda), v = k2 * sdLambda * sin(adLambda);
  areaRingSum.add(atan2(v, u));
  lambda0 = lambda, cosPhi0 = cosPhi, sinPhi0 = sinPhi;
}

// node_modules/d3-geo/src/cartesian.js
function spherical(cartesian3) {
  return [atan2(cartesian3[1], cartesian3[0]), asin(cartesian3[2])];
}
function cartesian(spherical3) {
  var lambda = spherical3[0], phi = spherical3[1], cosPhi = cos(phi);
  return [cosPhi * cos(lambda), cosPhi * sin(lambda), sin(phi)];
}
function cartesianDot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
function cartesianCross(a, b) {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}
function cartesianAddInPlace(a, b) {
  a[0] += b[0], a[1] += b[1], a[2] += b[2];
}
function cartesianScale(vector, k2) {
  return [vector[0] * k2, vector[1] * k2, vector[2] * k2];
}
function cartesianNormalizeInPlace(d) {
  var l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  d[0] /= l, d[1] /= l, d[2] /= l;
}

// node_modules/d3-geo/src/bounds.js
var lambda02;
var phi0;
var lambda1;
var phi1;
var lambda2;
var lambda002;
var phi002;
var p0;
var deltaSum;
var ranges;
var range2;
var boundsStream = {
  point: boundsPoint,
  lineStart: boundsLineStart,
  lineEnd: boundsLineEnd,
  polygonStart: function() {
    boundsStream.point = boundsRingPoint;
    boundsStream.lineStart = boundsRingStart;
    boundsStream.lineEnd = boundsRingEnd;
    deltaSum = new Adder();
    areaStream.polygonStart();
  },
  polygonEnd: function() {
    areaStream.polygonEnd();
    boundsStream.point = boundsPoint;
    boundsStream.lineStart = boundsLineStart;
    boundsStream.lineEnd = boundsLineEnd;
    if (areaRingSum < 0)
      lambda02 = -(lambda1 = 180), phi0 = -(phi1 = 90);
    else if (deltaSum > epsilon)
      phi1 = 90;
    else if (deltaSum < -epsilon)
      phi0 = -90;
    range2[0] = lambda02, range2[1] = lambda1;
  },
  sphere: function() {
    lambda02 = -(lambda1 = 180), phi0 = -(phi1 = 90);
  }
};
function boundsPoint(lambda, phi) {
  ranges.push(range2 = [lambda02 = lambda, lambda1 = lambda]);
  if (phi < phi0)
    phi0 = phi;
  if (phi > phi1)
    phi1 = phi;
}
function linePoint(lambda, phi) {
  var p = cartesian([lambda * radians, phi * radians]);
  if (p0) {
    var normal = cartesianCross(p0, p), equatorial = [normal[1], -normal[0], 0], inflection = cartesianCross(equatorial, normal);
    cartesianNormalizeInPlace(inflection);
    inflection = spherical(inflection);
    var delta = lambda - lambda2, sign3 = delta > 0 ? 1 : -1, lambdai = inflection[0] * degrees * sign3, phii, antimeridian = abs(delta) > 180;
    if (antimeridian ^ (sign3 * lambda2 < lambdai && lambdai < sign3 * lambda)) {
      phii = inflection[1] * degrees;
      if (phii > phi1)
        phi1 = phii;
    } else if (lambdai = (lambdai + 360) % 360 - 180, antimeridian ^ (sign3 * lambda2 < lambdai && lambdai < sign3 * lambda)) {
      phii = -inflection[1] * degrees;
      if (phii < phi0)
        phi0 = phii;
    } else {
      if (phi < phi0)
        phi0 = phi;
      if (phi > phi1)
        phi1 = phi;
    }
    if (antimeridian) {
      if (lambda < lambda2) {
        if (angle(lambda02, lambda) > angle(lambda02, lambda1))
          lambda1 = lambda;
      } else {
        if (angle(lambda, lambda1) > angle(lambda02, lambda1))
          lambda02 = lambda;
      }
    } else {
      if (lambda1 >= lambda02) {
        if (lambda < lambda02)
          lambda02 = lambda;
        if (lambda > lambda1)
          lambda1 = lambda;
      } else {
        if (lambda > lambda2) {
          if (angle(lambda02, lambda) > angle(lambda02, lambda1))
            lambda1 = lambda;
        } else {
          if (angle(lambda, lambda1) > angle(lambda02, lambda1))
            lambda02 = lambda;
        }
      }
    }
  } else {
    ranges.push(range2 = [lambda02 = lambda, lambda1 = lambda]);
  }
  if (phi < phi0)
    phi0 = phi;
  if (phi > phi1)
    phi1 = phi;
  p0 = p, lambda2 = lambda;
}
function boundsLineStart() {
  boundsStream.point = linePoint;
}
function boundsLineEnd() {
  range2[0] = lambda02, range2[1] = lambda1;
  boundsStream.point = boundsPoint;
  p0 = null;
}
function boundsRingPoint(lambda, phi) {
  if (p0) {
    var delta = lambda - lambda2;
    deltaSum.add(abs(delta) > 180 ? delta + (delta > 0 ? 360 : -360) : delta);
  } else {
    lambda002 = lambda, phi002 = phi;
  }
  areaStream.point(lambda, phi);
  linePoint(lambda, phi);
}
function boundsRingStart() {
  areaStream.lineStart();
}
function boundsRingEnd() {
  boundsRingPoint(lambda002, phi002);
  areaStream.lineEnd();
  if (abs(deltaSum) > epsilon)
    lambda02 = -(lambda1 = 180);
  range2[0] = lambda02, range2[1] = lambda1;
  p0 = null;
}
function angle(lambda03, lambda12) {
  return (lambda12 -= lambda03) < 0 ? lambda12 + 360 : lambda12;
}
function rangeCompare(a, b) {
  return a[0] - b[0];
}
function rangeContains(range3, x) {
  return range3[0] <= range3[1] ? range3[0] <= x && x <= range3[1] : x < range3[0] || range3[1] < x;
}
function bounds_default(feature) {
  var i, n, a, b, merged, deltaMax, delta;
  phi1 = lambda1 = -(lambda02 = phi0 = Infinity);
  ranges = [];
  stream_default(feature, boundsStream);
  if (n = ranges.length) {
    ranges.sort(rangeCompare);
    for (i = 1, a = ranges[0], merged = [a]; i < n; ++i) {
      b = ranges[i];
      if (rangeContains(a, b[0]) || rangeContains(a, b[1])) {
        if (angle(a[0], b[1]) > angle(a[0], a[1]))
          a[1] = b[1];
        if (angle(b[0], a[1]) > angle(a[0], a[1]))
          a[0] = b[0];
      } else {
        merged.push(a = b);
      }
    }
    for (deltaMax = -Infinity, n = merged.length - 1, i = 0, a = merged[n]; i <= n; a = b, ++i) {
      b = merged[i];
      if ((delta = angle(a[1], b[0])) > deltaMax)
        deltaMax = delta, lambda02 = b[0], lambda1 = a[1];
    }
  }
  ranges = range2 = null;
  return lambda02 === Infinity || phi0 === Infinity ? [[NaN, NaN], [NaN, NaN]] : [[lambda02, phi0], [lambda1, phi1]];
}

// node_modules/d3-geo/src/centroid.js
var W0;
var W1;
var X0;
var Y0;
var Z0;
var X1;
var Y1;
var Z1;
var X2;
var Y2;
var Z2;
var lambda003;
var phi003;
var x0;
var y0;
var z0;
var centroidStream = {
  sphere: noop,
  point: centroidPoint,
  lineStart: centroidLineStart,
  lineEnd: centroidLineEnd,
  polygonStart: function() {
    centroidStream.lineStart = centroidRingStart;
    centroidStream.lineEnd = centroidRingEnd;
  },
  polygonEnd: function() {
    centroidStream.lineStart = centroidLineStart;
    centroidStream.lineEnd = centroidLineEnd;
  }
};
function centroidPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi);
  centroidPointCartesian(cosPhi * cos(lambda), cosPhi * sin(lambda), sin(phi));
}
function centroidPointCartesian(x, y, z) {
  ++W0;
  X0 += (x - X0) / W0;
  Y0 += (y - Y0) / W0;
  Z0 += (z - Z0) / W0;
}
function centroidLineStart() {
  centroidStream.point = centroidLinePointFirst;
}
function centroidLinePointFirst(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi);
  x0 = cosPhi * cos(lambda);
  y0 = cosPhi * sin(lambda);
  z0 = sin(phi);
  centroidStream.point = centroidLinePoint;
  centroidPointCartesian(x0, y0, z0);
}
function centroidLinePoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi), x = cosPhi * cos(lambda), y = cosPhi * sin(lambda), z = sin(phi), w2 = atan2(sqrt((w2 = y0 * z - z0 * y) * w2 + (w2 = z0 * x - x0 * z) * w2 + (w2 = x0 * y - y0 * x) * w2), x0 * x + y0 * y + z0 * z);
  W1 += w2;
  X1 += w2 * (x0 + (x0 = x));
  Y1 += w2 * (y0 + (y0 = y));
  Z1 += w2 * (z0 + (z0 = z));
  centroidPointCartesian(x0, y0, z0);
}
function centroidLineEnd() {
  centroidStream.point = centroidPoint;
}
function centroidRingStart() {
  centroidStream.point = centroidRingPointFirst;
}
function centroidRingEnd() {
  centroidRingPoint(lambda003, phi003);
  centroidStream.point = centroidPoint;
}
function centroidRingPointFirst(lambda, phi) {
  lambda003 = lambda, phi003 = phi;
  lambda *= radians, phi *= radians;
  centroidStream.point = centroidRingPoint;
  var cosPhi = cos(phi);
  x0 = cosPhi * cos(lambda);
  y0 = cosPhi * sin(lambda);
  z0 = sin(phi);
  centroidPointCartesian(x0, y0, z0);
}
function centroidRingPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi), x = cosPhi * cos(lambda), y = cosPhi * sin(lambda), z = sin(phi), cx = y0 * z - z0 * y, cy = z0 * x - x0 * z, cz = x0 * y - y0 * x, m = hypot(cx, cy, cz), w2 = asin(m), v = m && -w2 / m;
  X2.add(v * cx);
  Y2.add(v * cy);
  Z2.add(v * cz);
  W1 += w2;
  X1 += w2 * (x0 + (x0 = x));
  Y1 += w2 * (y0 + (y0 = y));
  Z1 += w2 * (z0 + (z0 = z));
  centroidPointCartesian(x0, y0, z0);
}
function centroid_default(object) {
  W0 = W1 = X0 = Y0 = Z0 = X1 = Y1 = Z1 = 0;
  X2 = new Adder();
  Y2 = new Adder();
  Z2 = new Adder();
  stream_default(object, centroidStream);
  var x = +X2, y = +Y2, z = +Z2, m = hypot(x, y, z);
  if (m < epsilon2) {
    x = X1, y = Y1, z = Z1;
    if (W1 < epsilon)
      x = X0, y = Y0, z = Z0;
    m = hypot(x, y, z);
    if (m < epsilon2)
      return [NaN, NaN];
  }
  return [atan2(y, x) * degrees, asin(z / m) * degrees];
}

// node_modules/d3-geo/src/constant.js
function constant_default(x) {
  return function() {
    return x;
  };
}

// node_modules/d3-geo/src/compose.js
function compose_default(a, b) {
  function compose(x, y) {
    return x = a(x, y), b(x[0], x[1]);
  }
  if (a.invert && b.invert)
    compose.invert = function(x, y) {
      return x = b.invert(x, y), x && a.invert(x[0], x[1]);
    };
  return compose;
}

// node_modules/d3-geo/src/rotation.js
function rotationIdentity(lambda, phi) {
  if (abs(lambda) > pi)
    lambda -= Math.round(lambda / tau) * tau;
  return [lambda, phi];
}
rotationIdentity.invert = rotationIdentity;
function rotateRadians(deltaLambda, deltaPhi, deltaGamma) {
  return (deltaLambda %= tau) ? deltaPhi || deltaGamma ? compose_default(rotationLambda(deltaLambda), rotationPhiGamma(deltaPhi, deltaGamma)) : rotationLambda(deltaLambda) : deltaPhi || deltaGamma ? rotationPhiGamma(deltaPhi, deltaGamma) : rotationIdentity;
}
function forwardRotationLambda(deltaLambda) {
  return function(lambda, phi) {
    lambda += deltaLambda;
    if (abs(lambda) > pi)
      lambda -= Math.round(lambda / tau) * tau;
    return [lambda, phi];
  };
}
function rotationLambda(deltaLambda) {
  var rotation = forwardRotationLambda(deltaLambda);
  rotation.invert = forwardRotationLambda(-deltaLambda);
  return rotation;
}
function rotationPhiGamma(deltaPhi, deltaGamma) {
  var cosDeltaPhi = cos(deltaPhi), sinDeltaPhi = sin(deltaPhi), cosDeltaGamma = cos(deltaGamma), sinDeltaGamma = sin(deltaGamma);
  function rotation(lambda, phi) {
    var cosPhi = cos(phi), x = cos(lambda) * cosPhi, y = sin(lambda) * cosPhi, z = sin(phi), k2 = z * cosDeltaPhi + x * sinDeltaPhi;
    return [
      atan2(y * cosDeltaGamma - k2 * sinDeltaGamma, x * cosDeltaPhi - z * sinDeltaPhi),
      asin(k2 * cosDeltaGamma + y * sinDeltaGamma)
    ];
  }
  rotation.invert = function(lambda, phi) {
    var cosPhi = cos(phi), x = cos(lambda) * cosPhi, y = sin(lambda) * cosPhi, z = sin(phi), k2 = z * cosDeltaGamma - y * sinDeltaGamma;
    return [
      atan2(y * cosDeltaGamma + z * sinDeltaGamma, x * cosDeltaPhi + k2 * sinDeltaPhi),
      asin(k2 * cosDeltaPhi - x * sinDeltaPhi)
    ];
  };
  return rotation;
}
function rotation_default(rotate) {
  rotate = rotateRadians(rotate[0] * radians, rotate[1] * radians, rotate.length > 2 ? rotate[2] * radians : 0);
  function forward(coordinates) {
    coordinates = rotate(coordinates[0] * radians, coordinates[1] * radians);
    return coordinates[0] *= degrees, coordinates[1] *= degrees, coordinates;
  }
  forward.invert = function(coordinates) {
    coordinates = rotate.invert(coordinates[0] * radians, coordinates[1] * radians);
    return coordinates[0] *= degrees, coordinates[1] *= degrees, coordinates;
  };
  return forward;
}

// node_modules/d3-geo/src/circle.js
function circleStream(stream, radius, delta, direction, t0, t1) {
  if (!delta)
    return;
  var cosRadius = cos(radius), sinRadius = sin(radius), step = direction * delta;
  if (t0 == null) {
    t0 = radius + direction * tau;
    t1 = radius - step / 2;
  } else {
    t0 = circleRadius(cosRadius, t0);
    t1 = circleRadius(cosRadius, t1);
    if (direction > 0 ? t0 < t1 : t0 > t1)
      t0 += direction * tau;
  }
  for (var point, t = t0; direction > 0 ? t > t1 : t < t1; t -= step) {
    point = spherical([cosRadius, -sinRadius * cos(t), -sinRadius * sin(t)]);
    stream.point(point[0], point[1]);
  }
}
function circleRadius(cosRadius, point) {
  point = cartesian(point), point[0] -= cosRadius;
  cartesianNormalizeInPlace(point);
  var radius = acos(-point[1]);
  return ((-point[2] < 0 ? -radius : radius) + tau - epsilon) % tau;
}
function circle_default() {
  var center = constant_default([0, 0]), radius = constant_default(90), precision = constant_default(6), ring, rotate, stream = { point };
  function point(x, y) {
    ring.push(x = rotate(x, y));
    x[0] *= degrees, x[1] *= degrees;
  }
  function circle() {
    var c = center.apply(this, arguments), r = radius.apply(this, arguments) * radians, p = precision.apply(this, arguments) * radians;
    ring = [];
    rotate = rotateRadians(-c[0] * radians, -c[1] * radians, 0).invert;
    circleStream(stream, r, p, 1);
    c = { type: "Polygon", coordinates: [ring] };
    ring = rotate = null;
    return c;
  }
  circle.center = function(_) {
    return arguments.length ? (center = typeof _ === "function" ? _ : constant_default([+_[0], +_[1]]), circle) : center;
  };
  circle.radius = function(_) {
    return arguments.length ? (radius = typeof _ === "function" ? _ : constant_default(+_), circle) : radius;
  };
  circle.precision = function(_) {
    return arguments.length ? (precision = typeof _ === "function" ? _ : constant_default(+_), circle) : precision;
  };
  return circle;
}

// node_modules/d3-geo/src/clip/buffer.js
function buffer_default() {
  var lines2 = [], line;
  return {
    point: function(x, y, m) {
      line.push([x, y, m]);
    },
    lineStart: function() {
      lines2.push(line = []);
    },
    lineEnd: noop,
    rejoin: function() {
      if (lines2.length > 1)
        lines2.push(lines2.pop().concat(lines2.shift()));
    },
    result: function() {
      var result = lines2;
      lines2 = [];
      line = null;
      return result;
    }
  };
}

// node_modules/d3-geo/src/pointEqual.js
function pointEqual_default(a, b) {
  return abs(a[0] - b[0]) < epsilon && abs(a[1] - b[1]) < epsilon;
}

// node_modules/d3-geo/src/clip/rejoin.js
function Intersection(point, points2, other, entry) {
  this.x = point;
  this.z = points2;
  this.o = other;
  this.e = entry;
  this.v = false;
  this.n = this.p = null;
}
function rejoin_default(segments, compareIntersection2, startInside, interpolate, stream) {
  var subject = [], clip = [], i, n;
  segments.forEach(function(segment) {
    if ((n2 = segment.length - 1) <= 0)
      return;
    var n2, p02 = segment[0], p1 = segment[n2], x;
    if (pointEqual_default(p02, p1)) {
      if (!p02[2] && !p1[2]) {
        stream.lineStart();
        for (i = 0; i < n2; ++i)
          stream.point((p02 = segment[i])[0], p02[1]);
        stream.lineEnd();
        return;
      }
      p1[0] += 2 * epsilon;
    }
    subject.push(x = new Intersection(p02, segment, null, true));
    clip.push(x.o = new Intersection(p02, null, x, false));
    subject.push(x = new Intersection(p1, segment, null, false));
    clip.push(x.o = new Intersection(p1, null, x, true));
  });
  if (!subject.length)
    return;
  clip.sort(compareIntersection2);
  link(subject);
  link(clip);
  for (i = 0, n = clip.length; i < n; ++i) {
    clip[i].e = startInside = !startInside;
  }
  var start = subject[0], points2, point;
  while (1) {
    var current = start, isSubject = true;
    while (current.v)
      if ((current = current.n) === start)
        return;
    points2 = current.z;
    stream.lineStart();
    do {
      current.v = current.o.v = true;
      if (current.e) {
        if (isSubject) {
          for (i = 0, n = points2.length; i < n; ++i)
            stream.point((point = points2[i])[0], point[1]);
        } else {
          interpolate(current.x, current.n.x, 1, stream);
        }
        current = current.n;
      } else {
        if (isSubject) {
          points2 = current.p.z;
          for (i = points2.length - 1; i >= 0; --i)
            stream.point((point = points2[i])[0], point[1]);
        } else {
          interpolate(current.x, current.p.x, -1, stream);
        }
        current = current.p;
      }
      current = current.o;
      points2 = current.z;
      isSubject = !isSubject;
    } while (!current.v);
    stream.lineEnd();
  }
}
function link(array) {
  if (!(n = array.length))
    return;
  var n, i = 0, a = array[0], b;
  while (++i < n) {
    a.n = b = array[i];
    b.p = a;
    a = b;
  }
  a.n = b = array[0];
  b.p = a;
}

// node_modules/d3-geo/src/polygonContains.js
function longitude(point) {
  return abs(point[0]) <= pi ? point[0] : sign(point[0]) * ((abs(point[0]) + pi) % tau - pi);
}
function polygonContains_default(polygon, point) {
  var lambda = longitude(point), phi = point[1], sinPhi = sin(phi), normal = [sin(lambda), -cos(lambda), 0], angle3 = 0, winding = 0;
  var sum = new Adder();
  if (sinPhi === 1)
    phi = halfPi + epsilon;
  else if (sinPhi === -1)
    phi = -halfPi - epsilon;
  for (var i = 0, n = polygon.length; i < n; ++i) {
    if (!(m = (ring = polygon[i]).length))
      continue;
    var ring, m, point0 = ring[m - 1], lambda03 = longitude(point0), phi03 = point0[1] / 2 + quarterPi, sinPhi02 = sin(phi03), cosPhi02 = cos(phi03);
    for (var j = 0; j < m; ++j, lambda03 = lambda12, sinPhi02 = sinPhi1, cosPhi02 = cosPhi1, point0 = point1) {
      var point1 = ring[j], lambda12 = longitude(point1), phi12 = point1[1] / 2 + quarterPi, sinPhi1 = sin(phi12), cosPhi1 = cos(phi12), delta = lambda12 - lambda03, sign3 = delta >= 0 ? 1 : -1, absDelta = sign3 * delta, antimeridian = absDelta > pi, k2 = sinPhi02 * sinPhi1;
      sum.add(atan2(k2 * sign3 * sin(absDelta), cosPhi02 * cosPhi1 + k2 * cos(absDelta)));
      angle3 += antimeridian ? delta + sign3 * tau : delta;
      if (antimeridian ^ lambda03 >= lambda ^ lambda12 >= lambda) {
        var arc = cartesianCross(cartesian(point0), cartesian(point1));
        cartesianNormalizeInPlace(arc);
        var intersection = cartesianCross(normal, arc);
        cartesianNormalizeInPlace(intersection);
        var phiArc = (antimeridian ^ delta >= 0 ? -1 : 1) * asin(intersection[2]);
        if (phi > phiArc || phi === phiArc && (arc[0] || arc[1])) {
          winding += antimeridian ^ delta >= 0 ? 1 : -1;
        }
      }
    }
  }
  return (angle3 < -epsilon || angle3 < epsilon && sum < -epsilon2) ^ winding & 1;
}

// node_modules/d3-geo/src/clip/index.js
function clip_default(pointVisible, clipLine, interpolate, start) {
  return function(sink) {
    var line = clipLine(sink), ringBuffer = buffer_default(), ringSink = clipLine(ringBuffer), polygonStarted = false, polygon, segments, ring;
    var clip = {
      point,
      lineStart,
      lineEnd,
      polygonStart: function() {
        clip.point = pointRing;
        clip.lineStart = ringStart;
        clip.lineEnd = ringEnd;
        segments = [];
        polygon = [];
      },
      polygonEnd: function() {
        clip.point = point;
        clip.lineStart = lineStart;
        clip.lineEnd = lineEnd;
        segments = merge(segments);
        var startInside = polygonContains_default(polygon, start);
        if (segments.length) {
          if (!polygonStarted)
            sink.polygonStart(), polygonStarted = true;
          rejoin_default(segments, compareIntersection, startInside, interpolate, sink);
        } else if (startInside) {
          if (!polygonStarted)
            sink.polygonStart(), polygonStarted = true;
          sink.lineStart();
          interpolate(null, null, 1, sink);
          sink.lineEnd();
        }
        if (polygonStarted)
          sink.polygonEnd(), polygonStarted = false;
        segments = polygon = null;
      },
      sphere: function() {
        sink.polygonStart();
        sink.lineStart();
        interpolate(null, null, 1, sink);
        sink.lineEnd();
        sink.polygonEnd();
      }
    };
    function point(lambda, phi) {
      if (pointVisible(lambda, phi))
        sink.point(lambda, phi);
    }
    function pointLine(lambda, phi) {
      line.point(lambda, phi);
    }
    function lineStart() {
      clip.point = pointLine;
      line.lineStart();
    }
    function lineEnd() {
      clip.point = point;
      line.lineEnd();
    }
    function pointRing(lambda, phi) {
      ring.push([lambda, phi]);
      ringSink.point(lambda, phi);
    }
    function ringStart() {
      ringSink.lineStart();
      ring = [];
    }
    function ringEnd() {
      pointRing(ring[0][0], ring[0][1]);
      ringSink.lineEnd();
      var clean = ringSink.clean(), ringSegments = ringBuffer.result(), i, n = ringSegments.length, m, segment, point2;
      ring.pop();
      polygon.push(ring);
      ring = null;
      if (!n)
        return;
      if (clean & 1) {
        segment = ringSegments[0];
        if ((m = segment.length - 1) > 0) {
          if (!polygonStarted)
            sink.polygonStart(), polygonStarted = true;
          sink.lineStart();
          for (i = 0; i < m; ++i)
            sink.point((point2 = segment[i])[0], point2[1]);
          sink.lineEnd();
        }
        return;
      }
      if (n > 1 && clean & 2)
        ringSegments.push(ringSegments.pop().concat(ringSegments.shift()));
      segments.push(ringSegments.filter(validSegment));
    }
    return clip;
  };
}
function validSegment(segment) {
  return segment.length > 1;
}
function compareIntersection(a, b) {
  return ((a = a.x)[0] < 0 ? a[1] - halfPi - epsilon : halfPi - a[1]) - ((b = b.x)[0] < 0 ? b[1] - halfPi - epsilon : halfPi - b[1]);
}

// node_modules/d3-geo/src/clip/antimeridian.js
var antimeridian_default = clip_default(
  function() {
    return true;
  },
  clipAntimeridianLine,
  clipAntimeridianInterpolate,
  [-pi, -halfPi]
);
function clipAntimeridianLine(stream) {
  var lambda03 = NaN, phi03 = NaN, sign0 = NaN, clean;
  return {
    lineStart: function() {
      stream.lineStart();
      clean = 1;
    },
    point: function(lambda12, phi12) {
      var sign1 = lambda12 > 0 ? pi : -pi, delta = abs(lambda12 - lambda03);
      if (abs(delta - pi) < epsilon) {
        stream.point(lambda03, phi03 = (phi03 + phi12) / 2 > 0 ? halfPi : -halfPi);
        stream.point(sign0, phi03);
        stream.lineEnd();
        stream.lineStart();
        stream.point(sign1, phi03);
        stream.point(lambda12, phi03);
        clean = 0;
      } else if (sign0 !== sign1 && delta >= pi) {
        if (abs(lambda03 - sign0) < epsilon)
          lambda03 -= sign0 * epsilon;
        if (abs(lambda12 - sign1) < epsilon)
          lambda12 -= sign1 * epsilon;
        phi03 = clipAntimeridianIntersect(lambda03, phi03, lambda12, phi12);
        stream.point(sign0, phi03);
        stream.lineEnd();
        stream.lineStart();
        stream.point(sign1, phi03);
        clean = 0;
      }
      stream.point(lambda03 = lambda12, phi03 = phi12);
      sign0 = sign1;
    },
    lineEnd: function() {
      stream.lineEnd();
      lambda03 = phi03 = NaN;
    },
    clean: function() {
      return 2 - clean;
    }
  };
}
function clipAntimeridianIntersect(lambda03, phi03, lambda12, phi12) {
  var cosPhi02, cosPhi1, sinLambda0Lambda1 = sin(lambda03 - lambda12);
  return abs(sinLambda0Lambda1) > epsilon ? atan((sin(phi03) * (cosPhi1 = cos(phi12)) * sin(lambda12) - sin(phi12) * (cosPhi02 = cos(phi03)) * sin(lambda03)) / (cosPhi02 * cosPhi1 * sinLambda0Lambda1)) : (phi03 + phi12) / 2;
}
function clipAntimeridianInterpolate(from, to, direction, stream) {
  var phi;
  if (from == null) {
    phi = direction * halfPi;
    stream.point(-pi, phi);
    stream.point(0, phi);
    stream.point(pi, phi);
    stream.point(pi, 0);
    stream.point(pi, -phi);
    stream.point(0, -phi);
    stream.point(-pi, -phi);
    stream.point(-pi, 0);
    stream.point(-pi, phi);
  } else if (abs(from[0] - to[0]) > epsilon) {
    var lambda = from[0] < to[0] ? pi : -pi;
    phi = direction * lambda / 2;
    stream.point(-lambda, phi);
    stream.point(0, phi);
    stream.point(lambda, phi);
  } else {
    stream.point(to[0], to[1]);
  }
}

// node_modules/d3-geo/src/clip/circle.js
function circle_default2(radius) {
  var cr = cos(radius), delta = 6 * radians, smallRadius = cr > 0, notHemisphere = abs(cr) > epsilon;
  function interpolate(from, to, direction, stream) {
    circleStream(stream, radius, delta, direction, from, to);
  }
  function visible(lambda, phi) {
    return cos(lambda) * cos(phi) > cr;
  }
  function clipLine(stream) {
    var point0, c0, v0, v00, clean;
    return {
      lineStart: function() {
        v00 = v0 = false;
        clean = 1;
      },
      point: function(lambda, phi) {
        var point1 = [lambda, phi], point2, v = visible(lambda, phi), c = smallRadius ? v ? 0 : code(lambda, phi) : v ? code(lambda + (lambda < 0 ? pi : -pi), phi) : 0;
        if (!point0 && (v00 = v0 = v))
          stream.lineStart();
        if (v !== v0) {
          point2 = intersect(point0, point1);
          if (!point2 || pointEqual_default(point0, point2) || pointEqual_default(point1, point2))
            point1[2] = 1;
        }
        if (v !== v0) {
          clean = 0;
          if (v) {
            stream.lineStart();
            point2 = intersect(point1, point0);
            stream.point(point2[0], point2[1]);
          } else {
            point2 = intersect(point0, point1);
            stream.point(point2[0], point2[1], 2);
            stream.lineEnd();
          }
          point0 = point2;
        } else if (notHemisphere && point0 && smallRadius ^ v) {
          var t;
          if (!(c & c0) && (t = intersect(point1, point0, true))) {
            clean = 0;
            if (smallRadius) {
              stream.lineStart();
              stream.point(t[0][0], t[0][1]);
              stream.point(t[1][0], t[1][1]);
              stream.lineEnd();
            } else {
              stream.point(t[1][0], t[1][1]);
              stream.lineEnd();
              stream.lineStart();
              stream.point(t[0][0], t[0][1], 3);
            }
          }
        }
        if (v && (!point0 || !pointEqual_default(point0, point1))) {
          stream.point(point1[0], point1[1]);
        }
        point0 = point1, v0 = v, c0 = c;
      },
      lineEnd: function() {
        if (v0)
          stream.lineEnd();
        point0 = null;
      },
      // Rejoin first and last segments if there were intersections and the first
      // and last points were visible.
      clean: function() {
        return clean | (v00 && v0) << 1;
      }
    };
  }
  function intersect(a, b, two) {
    var pa = cartesian(a), pb = cartesian(b);
    var n1 = [1, 0, 0], n2 = cartesianCross(pa, pb), n2n2 = cartesianDot(n2, n2), n1n2 = n2[0], determinant = n2n2 - n1n2 * n1n2;
    if (!determinant)
      return !two && a;
    var c1 = cr * n2n2 / determinant, c2 = -cr * n1n2 / determinant, n1xn2 = cartesianCross(n1, n2), A5 = cartesianScale(n1, c1), B2 = cartesianScale(n2, c2);
    cartesianAddInPlace(A5, B2);
    var u = n1xn2, w2 = cartesianDot(A5, u), uu = cartesianDot(u, u), t2 = w2 * w2 - uu * (cartesianDot(A5, A5) - 1);
    if (t2 < 0)
      return;
    var t = sqrt(t2), q = cartesianScale(u, (-w2 - t) / uu);
    cartesianAddInPlace(q, A5);
    q = spherical(q);
    if (!two)
      return q;
    var lambda03 = a[0], lambda12 = b[0], phi03 = a[1], phi12 = b[1], z;
    if (lambda12 < lambda03)
      z = lambda03, lambda03 = lambda12, lambda12 = z;
    var delta2 = lambda12 - lambda03, polar = abs(delta2 - pi) < epsilon, meridian = polar || delta2 < epsilon;
    if (!polar && phi12 < phi03)
      z = phi03, phi03 = phi12, phi12 = z;
    if (meridian ? polar ? phi03 + phi12 > 0 ^ q[1] < (abs(q[0] - lambda03) < epsilon ? phi03 : phi12) : phi03 <= q[1] && q[1] <= phi12 : delta2 > pi ^ (lambda03 <= q[0] && q[0] <= lambda12)) {
      var q1 = cartesianScale(u, (-w2 + t) / uu);
      cartesianAddInPlace(q1, A5);
      return [q, spherical(q1)];
    }
  }
  function code(lambda, phi) {
    var r = smallRadius ? radius : pi - radius, code2 = 0;
    if (lambda < -r)
      code2 |= 1;
    else if (lambda > r)
      code2 |= 2;
    if (phi < -r)
      code2 |= 4;
    else if (phi > r)
      code2 |= 8;
    return code2;
  }
  return clip_default(visible, clipLine, interpolate, smallRadius ? [0, -radius] : [-pi, radius - pi]);
}

// node_modules/d3-geo/src/clip/line.js
function line_default(a, b, x03, y03, x12, y12) {
  var ax = a[0], ay = a[1], bx = b[0], by = b[1], t0 = 0, t1 = 1, dx = bx - ax, dy = by - ay, r;
  r = x03 - ax;
  if (!dx && r > 0)
    return;
  r /= dx;
  if (dx < 0) {
    if (r < t0)
      return;
    if (r < t1)
      t1 = r;
  } else if (dx > 0) {
    if (r > t1)
      return;
    if (r > t0)
      t0 = r;
  }
  r = x12 - ax;
  if (!dx && r < 0)
    return;
  r /= dx;
  if (dx < 0) {
    if (r > t1)
      return;
    if (r > t0)
      t0 = r;
  } else if (dx > 0) {
    if (r < t0)
      return;
    if (r < t1)
      t1 = r;
  }
  r = y03 - ay;
  if (!dy && r > 0)
    return;
  r /= dy;
  if (dy < 0) {
    if (r < t0)
      return;
    if (r < t1)
      t1 = r;
  } else if (dy > 0) {
    if (r > t1)
      return;
    if (r > t0)
      t0 = r;
  }
  r = y12 - ay;
  if (!dy && r < 0)
    return;
  r /= dy;
  if (dy < 0) {
    if (r > t1)
      return;
    if (r > t0)
      t0 = r;
  } else if (dy > 0) {
    if (r < t0)
      return;
    if (r < t1)
      t1 = r;
  }
  if (t0 > 0)
    a[0] = ax + t0 * dx, a[1] = ay + t0 * dy;
  if (t1 < 1)
    b[0] = ax + t1 * dx, b[1] = ay + t1 * dy;
  return true;
}

// node_modules/d3-geo/src/clip/rectangle.js
var clipMax = 1e9;
var clipMin = -clipMax;
function clipRectangle(x03, y03, x12, y12) {
  function visible(x, y) {
    return x03 <= x && x <= x12 && y03 <= y && y <= y12;
  }
  function interpolate(from, to, direction, stream) {
    var a = 0, a1 = 0;
    if (from == null || (a = corner(from, direction)) !== (a1 = corner(to, direction)) || comparePoint(from, to) < 0 ^ direction > 0) {
      do
        stream.point(a === 0 || a === 3 ? x03 : x12, a > 1 ? y12 : y03);
      while ((a = (a + direction + 4) % 4) !== a1);
    } else {
      stream.point(to[0], to[1]);
    }
  }
  function corner(p, direction) {
    return abs(p[0] - x03) < epsilon ? direction > 0 ? 0 : 3 : abs(p[0] - x12) < epsilon ? direction > 0 ? 2 : 1 : abs(p[1] - y03) < epsilon ? direction > 0 ? 1 : 0 : direction > 0 ? 3 : 2;
  }
  function compareIntersection2(a, b) {
    return comparePoint(a.x, b.x);
  }
  function comparePoint(a, b) {
    var ca = corner(a, 1), cb = corner(b, 1);
    return ca !== cb ? ca - cb : ca === 0 ? b[1] - a[1] : ca === 1 ? a[0] - b[0] : ca === 2 ? a[1] - b[1] : b[0] - a[0];
  }
  return function(stream) {
    var activeStream = stream, bufferStream = buffer_default(), segments, polygon, ring, x__, y__, v__, x_, y_, v_, first, clean;
    var clipStream = {
      point,
      lineStart,
      lineEnd,
      polygonStart,
      polygonEnd
    };
    function point(x, y) {
      if (visible(x, y))
        activeStream.point(x, y);
    }
    function polygonInside() {
      var winding = 0;
      for (var i = 0, n = polygon.length; i < n; ++i) {
        for (var ring2 = polygon[i], j = 1, m = ring2.length, point2 = ring2[0], a0, a1, b0 = point2[0], b1 = point2[1]; j < m; ++j) {
          a0 = b0, a1 = b1, point2 = ring2[j], b0 = point2[0], b1 = point2[1];
          if (a1 <= y12) {
            if (b1 > y12 && (b0 - a0) * (y12 - a1) > (b1 - a1) * (x03 - a0))
              ++winding;
          } else {
            if (b1 <= y12 && (b0 - a0) * (y12 - a1) < (b1 - a1) * (x03 - a0))
              --winding;
          }
        }
      }
      return winding;
    }
    function polygonStart() {
      activeStream = bufferStream, segments = [], polygon = [], clean = true;
    }
    function polygonEnd() {
      var startInside = polygonInside(), cleanInside = clean && startInside, visible2 = (segments = merge(segments)).length;
      if (cleanInside || visible2) {
        stream.polygonStart();
        if (cleanInside) {
          stream.lineStart();
          interpolate(null, null, 1, stream);
          stream.lineEnd();
        }
        if (visible2) {
          rejoin_default(segments, compareIntersection2, startInside, interpolate, stream);
        }
        stream.polygonEnd();
      }
      activeStream = stream, segments = polygon = ring = null;
    }
    function lineStart() {
      clipStream.point = linePoint2;
      if (polygon)
        polygon.push(ring = []);
      first = true;
      v_ = false;
      x_ = y_ = NaN;
    }
    function lineEnd() {
      if (segments) {
        linePoint2(x__, y__);
        if (v__ && v_)
          bufferStream.rejoin();
        segments.push(bufferStream.result());
      }
      clipStream.point = point;
      if (v_)
        activeStream.lineEnd();
    }
    function linePoint2(x, y) {
      var v = visible(x, y);
      if (polygon)
        ring.push([x, y]);
      if (first) {
        x__ = x, y__ = y, v__ = v;
        first = false;
        if (v) {
          activeStream.lineStart();
          activeStream.point(x, y);
        }
      } else {
        if (v && v_)
          activeStream.point(x, y);
        else {
          var a = [x_ = Math.max(clipMin, Math.min(clipMax, x_)), y_ = Math.max(clipMin, Math.min(clipMax, y_))], b = [x = Math.max(clipMin, Math.min(clipMax, x)), y = Math.max(clipMin, Math.min(clipMax, y))];
          if (line_default(a, b, x03, y03, x12, y12)) {
            if (!v_) {
              activeStream.lineStart();
              activeStream.point(a[0], a[1]);
            }
            activeStream.point(b[0], b[1]);
            if (!v)
              activeStream.lineEnd();
            clean = false;
          } else if (v) {
            activeStream.lineStart();
            activeStream.point(x, y);
            clean = false;
          }
        }
      }
      x_ = x, y_ = y, v_ = v;
    }
    return clipStream;
  };
}

// node_modules/d3-geo/src/interpolate.js
function interpolate_default(a, b) {
  var x03 = a[0] * radians, y03 = a[1] * radians, x12 = b[0] * radians, y12 = b[1] * radians, cy0 = cos(y03), sy0 = sin(y03), cy1 = cos(y12), sy1 = sin(y12), kx0 = cy0 * cos(x03), ky0 = cy0 * sin(x03), kx1 = cy1 * cos(x12), ky1 = cy1 * sin(x12), d = 2 * asin(sqrt(haversin(y12 - y03) + cy0 * cy1 * haversin(x12 - x03))), k2 = sin(d);
  var interpolate = d ? function(t) {
    var B2 = sin(t *= d) / k2, A5 = sin(d - t) / k2, x = A5 * kx0 + B2 * kx1, y = A5 * ky0 + B2 * ky1, z = A5 * sy0 + B2 * sy1;
    return [
      atan2(y, x) * degrees,
      atan2(z, sqrt(x * x + y * y)) * degrees
    ];
  } : function() {
    return [x03 * degrees, y03 * degrees];
  };
  interpolate.distance = d;
  return interpolate;
}

// node_modules/d3-geo/src/identity.js
var identity_default = (x) => x;

// node_modules/d3-geo/src/path/bounds.js
var x02 = Infinity;
var y02 = x02;
var x1 = -x02;
var y1 = x1;
var boundsStream2 = {
  point: boundsPoint2,
  lineStart: noop,
  lineEnd: noop,
  polygonStart: noop,
  polygonEnd: noop,
  result: function() {
    var bounds = [[x02, y02], [x1, y1]];
    x1 = y1 = -(y02 = x02 = Infinity);
    return bounds;
  }
};
function boundsPoint2(x, y) {
  if (x < x02)
    x02 = x;
  if (x > x1)
    x1 = x;
  if (y < y02)
    y02 = y;
  if (y > y1)
    y1 = y;
}
var bounds_default2 = boundsStream2;

// node_modules/d3-geo/src/transform.js
function transform_default(methods) {
  return {
    stream: transformer(methods)
  };
}
function transformer(methods) {
  return function(stream) {
    var s = new TransformStream();
    for (var key in methods)
      s[key] = methods[key];
    s.stream = stream;
    return s;
  };
}
function TransformStream() {
}
TransformStream.prototype = {
  constructor: TransformStream,
  point: function(x, y) {
    this.stream.point(x, y);
  },
  sphere: function() {
    this.stream.sphere();
  },
  lineStart: function() {
    this.stream.lineStart();
  },
  lineEnd: function() {
    this.stream.lineEnd();
  },
  polygonStart: function() {
    this.stream.polygonStart();
  },
  polygonEnd: function() {
    this.stream.polygonEnd();
  }
};

// node_modules/d3-geo/src/projection/fit.js
function fit(projection2, fitBounds, object) {
  var clip = projection2.clipExtent && projection2.clipExtent();
  projection2.scale(150).translate([0, 0]);
  if (clip != null)
    projection2.clipExtent(null);
  stream_default(object, projection2.stream(bounds_default2));
  fitBounds(bounds_default2.result());
  if (clip != null)
    projection2.clipExtent(clip);
  return projection2;
}
function fitExtent(projection2, extent, object) {
  return fit(projection2, function(b) {
    var w2 = extent[1][0] - extent[0][0], h = extent[1][1] - extent[0][1], k2 = Math.min(w2 / (b[1][0] - b[0][0]), h / (b[1][1] - b[0][1])), x = +extent[0][0] + (w2 - k2 * (b[1][0] + b[0][0])) / 2, y = +extent[0][1] + (h - k2 * (b[1][1] + b[0][1])) / 2;
    projection2.scale(150 * k2).translate([x, y]);
  }, object);
}
function fitSize(projection2, size, object) {
  return fitExtent(projection2, [[0, 0], size], object);
}
function fitWidth(projection2, width, object) {
  return fit(projection2, function(b) {
    var w2 = +width, k2 = w2 / (b[1][0] - b[0][0]), x = (w2 - k2 * (b[1][0] + b[0][0])) / 2, y = -k2 * b[0][1];
    projection2.scale(150 * k2).translate([x, y]);
  }, object);
}
function fitHeight(projection2, height, object) {
  return fit(projection2, function(b) {
    var h = +height, k2 = h / (b[1][1] - b[0][1]), x = -k2 * b[0][0], y = (h - k2 * (b[1][1] + b[0][1])) / 2;
    projection2.scale(150 * k2).translate([x, y]);
  }, object);
}

// node_modules/d3-geo/src/projection/resample.js
var maxDepth = 16;
var cosMinDistance = cos(30 * radians);
function resample_default(project, delta2) {
  return +delta2 ? resample(project, delta2) : resampleNone(project);
}
function resampleNone(project) {
  return transformer({
    point: function(x, y) {
      x = project(x, y);
      this.stream.point(x[0], x[1]);
    }
  });
}
function resample(project, delta2) {
  function resampleLineTo(x03, y03, lambda03, a0, b0, c0, x12, y12, lambda12, a1, b1, c1, depth, stream) {
    var dx = x12 - x03, dy = y12 - y03, d2 = dx * dx + dy * dy;
    if (d2 > 4 * delta2 && depth--) {
      var a = a0 + a1, b = b0 + b1, c = c0 + c1, m = sqrt(a * a + b * b + c * c), phi2 = asin(c /= m), lambda22 = abs(abs(c) - 1) < epsilon || abs(lambda03 - lambda12) < epsilon ? (lambda03 + lambda12) / 2 : atan2(b, a), p = project(lambda22, phi2), x2 = p[0], y2 = p[1], dx2 = x2 - x03, dy2 = y2 - y03, dz = dy * dx2 - dx * dy2;
      if (dz * dz / d2 > delta2 || abs((dx * dx2 + dy * dy2) / d2 - 0.5) > 0.3 || a0 * a1 + b0 * b1 + c0 * c1 < cosMinDistance) {
        resampleLineTo(x03, y03, lambda03, a0, b0, c0, x2, y2, lambda22, a /= m, b /= m, c, depth, stream);
        stream.point(x2, y2);
        resampleLineTo(x2, y2, lambda22, a, b, c, x12, y12, lambda12, a1, b1, c1, depth, stream);
      }
    }
  }
  return function(stream) {
    var lambda004, x00, y00, a00, b00, c00, lambda03, x03, y03, a0, b0, c0;
    var resampleStream = {
      point,
      lineStart,
      lineEnd,
      polygonStart: function() {
        stream.polygonStart();
        resampleStream.lineStart = ringStart;
      },
      polygonEnd: function() {
        stream.polygonEnd();
        resampleStream.lineStart = lineStart;
      }
    };
    function point(x, y) {
      x = project(x, y);
      stream.point(x[0], x[1]);
    }
    function lineStart() {
      x03 = NaN;
      resampleStream.point = linePoint2;
      stream.lineStart();
    }
    function linePoint2(lambda, phi) {
      var c = cartesian([lambda, phi]), p = project(lambda, phi);
      resampleLineTo(x03, y03, lambda03, a0, b0, c0, x03 = p[0], y03 = p[1], lambda03 = lambda, a0 = c[0], b0 = c[1], c0 = c[2], maxDepth, stream);
      stream.point(x03, y03);
    }
    function lineEnd() {
      resampleStream.point = point;
      stream.lineEnd();
    }
    function ringStart() {
      lineStart();
      resampleStream.point = ringPoint;
      resampleStream.lineEnd = ringEnd;
    }
    function ringPoint(lambda, phi) {
      linePoint2(lambda004 = lambda, phi), x00 = x03, y00 = y03, a00 = a0, b00 = b0, c00 = c0;
      resampleStream.point = linePoint2;
    }
    function ringEnd() {
      resampleLineTo(x03, y03, lambda03, a0, b0, c0, x00, y00, lambda004, a00, b00, c00, maxDepth, stream);
      resampleStream.lineEnd = lineEnd;
      lineEnd();
    }
    return resampleStream;
  };
}

// node_modules/d3-geo/src/projection/index.js
var transformRadians = transformer({
  point: function(x, y) {
    this.stream.point(x * radians, y * radians);
  }
});
function transformRotate(rotate) {
  return transformer({
    point: function(x, y) {
      var r = rotate(x, y);
      return this.stream.point(r[0], r[1]);
    }
  });
}
function scaleTranslate(k2, dx, dy, sx, sy) {
  function transform(x, y) {
    x *= sx;
    y *= sy;
    return [dx + k2 * x, dy - k2 * y];
  }
  transform.invert = function(x, y) {
    return [(x - dx) / k2 * sx, (dy - y) / k2 * sy];
  };
  return transform;
}
function scaleTranslateRotate(k2, dx, dy, sx, sy, alpha) {
  if (!alpha)
    return scaleTranslate(k2, dx, dy, sx, sy);
  var cosAlpha = cos(alpha), sinAlpha = sin(alpha), a = cosAlpha * k2, b = sinAlpha * k2, ai = cosAlpha / k2, bi = sinAlpha / k2, ci = (sinAlpha * dy - cosAlpha * dx) / k2, fi = (sinAlpha * dx + cosAlpha * dy) / k2;
  function transform(x, y) {
    x *= sx;
    y *= sy;
    return [a * x - b * y + dx, dy - b * x - a * y];
  }
  transform.invert = function(x, y) {
    return [sx * (ai * x - bi * y + ci), sy * (fi - bi * x - ai * y)];
  };
  return transform;
}
function projection(project) {
  return projectionMutator(function() {
    return project;
  })();
}
function projectionMutator(projectAt) {
  var project, k2 = 150, x = 480, y = 250, lambda = 0, phi = 0, deltaLambda = 0, deltaPhi = 0, deltaGamma = 0, rotate, alpha = 0, sx = 1, sy = 1, theta = null, preclip = antimeridian_default, x03 = null, y03, x12, y12, postclip = identity_default, delta2 = 0.5, projectResample, projectTransform, projectRotateTransform, cache, cacheStream;
  function projection2(point) {
    return projectRotateTransform(point[0] * radians, point[1] * radians);
  }
  function invert(point) {
    point = projectRotateTransform.invert(point[0], point[1]);
    return point && [point[0] * degrees, point[1] * degrees];
  }
  projection2.stream = function(stream) {
    return cache && cacheStream === stream ? cache : cache = transformRadians(transformRotate(rotate)(preclip(projectResample(postclip(cacheStream = stream)))));
  };
  projection2.preclip = function(_) {
    return arguments.length ? (preclip = _, theta = void 0, reset()) : preclip;
  };
  projection2.postclip = function(_) {
    return arguments.length ? (postclip = _, x03 = y03 = x12 = y12 = null, reset()) : postclip;
  };
  projection2.clipAngle = function(_) {
    return arguments.length ? (preclip = +_ ? circle_default2(theta = _ * radians) : (theta = null, antimeridian_default), reset()) : theta * degrees;
  };
  projection2.clipExtent = function(_) {
    return arguments.length ? (postclip = _ == null ? (x03 = y03 = x12 = y12 = null, identity_default) : clipRectangle(x03 = +_[0][0], y03 = +_[0][1], x12 = +_[1][0], y12 = +_[1][1]), reset()) : x03 == null ? null : [[x03, y03], [x12, y12]];
  };
  projection2.scale = function(_) {
    return arguments.length ? (k2 = +_, recenter()) : k2;
  };
  projection2.translate = function(_) {
    return arguments.length ? (x = +_[0], y = +_[1], recenter()) : [x, y];
  };
  projection2.center = function(_) {
    return arguments.length ? (lambda = _[0] % 360 * radians, phi = _[1] % 360 * radians, recenter()) : [lambda * degrees, phi * degrees];
  };
  projection2.rotate = function(_) {
    return arguments.length ? (deltaLambda = _[0] % 360 * radians, deltaPhi = _[1] % 360 * radians, deltaGamma = _.length > 2 ? _[2] % 360 * radians : 0, recenter()) : [deltaLambda * degrees, deltaPhi * degrees, deltaGamma * degrees];
  };
  projection2.angle = function(_) {
    return arguments.length ? (alpha = _ % 360 * radians, recenter()) : alpha * degrees;
  };
  projection2.reflectX = function(_) {
    return arguments.length ? (sx = _ ? -1 : 1, recenter()) : sx < 0;
  };
  projection2.reflectY = function(_) {
    return arguments.length ? (sy = _ ? -1 : 1, recenter()) : sy < 0;
  };
  projection2.precision = function(_) {
    return arguments.length ? (projectResample = resample_default(projectTransform, delta2 = _ * _), reset()) : sqrt(delta2);
  };
  projection2.fitExtent = function(extent, object) {
    return fitExtent(projection2, extent, object);
  };
  projection2.fitSize = function(size, object) {
    return fitSize(projection2, size, object);
  };
  projection2.fitWidth = function(width, object) {
    return fitWidth(projection2, width, object);
  };
  projection2.fitHeight = function(height, object) {
    return fitHeight(projection2, height, object);
  };
  function recenter() {
    var center = scaleTranslateRotate(k2, 0, 0, sx, sy, alpha).apply(null, project(lambda, phi)), transform = scaleTranslateRotate(k2, x - center[0], y - center[1], sx, sy, alpha);
    rotate = rotateRadians(deltaLambda, deltaPhi, deltaGamma);
    projectTransform = compose_default(project, transform);
    projectRotateTransform = compose_default(rotate, projectTransform);
    projectResample = resample_default(projectTransform, delta2);
    return reset();
  }
  function reset() {
    cache = cacheStream = null;
    return projection2;
  }
  return function() {
    project = projectAt.apply(this, arguments);
    projection2.invert = project.invert && invert;
    return recenter();
  };
}

// node_modules/d3-geo/src/projection/conic.js
function conicProjection(projectAt) {
  var phi03 = 0, phi12 = pi / 3, m = projectionMutator(projectAt), p = m(phi03, phi12);
  p.parallels = function(_) {
    return arguments.length ? m(phi03 = _[0] * radians, phi12 = _[1] * radians) : [phi03 * degrees, phi12 * degrees];
  };
  return p;
}

// node_modules/d3-geo/src/projection/cylindricalEqualArea.js
function cylindricalEqualAreaRaw(phi03) {
  var cosPhi02 = cos(phi03);
  function forward(lambda, phi) {
    return [lambda * cosPhi02, sin(phi) / cosPhi02];
  }
  forward.invert = function(x, y) {
    return [x / cosPhi02, asin(y * cosPhi02)];
  };
  return forward;
}

// node_modules/d3-geo/src/projection/conicEqualArea.js
function conicEqualAreaRaw(y03, y12) {
  var sy0 = sin(y03), n = (sy0 + sin(y12)) / 2;
  if (abs(n) < epsilon)
    return cylindricalEqualAreaRaw(y03);
  var c = 1 + sy0 * (2 * n - sy0), r0 = sqrt(c) / n;
  function project(x, y) {
    var r = sqrt(c - 2 * n * sin(y)) / n;
    return [r * sin(x *= n), r0 - r * cos(x)];
  }
  project.invert = function(x, y) {
    var r0y = r0 - y, l = atan2(x, abs(r0y)) * sign(r0y);
    if (r0y * n < 0)
      l -= pi * sign(x) * sign(r0y);
    return [l / n, asin((c - (x * x + r0y * r0y) * n * n) / (2 * n))];
  };
  return project;
}
function conicEqualArea_default() {
  return conicProjection(conicEqualAreaRaw).scale(155.424).center([0, 33.6442]);
}

// node_modules/d3-geo/src/projection/albers.js
function albers_default() {
  return conicEqualArea_default().parallels([29.5, 45.5]).scale(1070).translate([480, 250]).rotate([96, 0]).center([-0.6, 38.7]);
}

// node_modules/d3-geo/src/projection/albersUsa.js
function multiplex(streams) {
  var n = streams.length;
  return {
    point: function(x, y) {
      var i = -1;
      while (++i < n)
        streams[i].point(x, y);
    },
    sphere: function() {
      var i = -1;
      while (++i < n)
        streams[i].sphere();
    },
    lineStart: function() {
      var i = -1;
      while (++i < n)
        streams[i].lineStart();
    },
    lineEnd: function() {
      var i = -1;
      while (++i < n)
        streams[i].lineEnd();
    },
    polygonStart: function() {
      var i = -1;
      while (++i < n)
        streams[i].polygonStart();
    },
    polygonEnd: function() {
      var i = -1;
      while (++i < n)
        streams[i].polygonEnd();
    }
  };
}
function albersUsa_default() {
  var cache, cacheStream, lower48 = albers_default(), lower48Point, alaska = conicEqualArea_default().rotate([154, 0]).center([-2, 58.5]).parallels([55, 65]), alaskaPoint, hawaii = conicEqualArea_default().rotate([157, 0]).center([-3, 19.9]).parallels([8, 18]), hawaiiPoint, point, pointStream = { point: function(x, y) {
    point = [x, y];
  } };
  function albersUsa(coordinates) {
    var x = coordinates[0], y = coordinates[1];
    return point = null, (lower48Point.point(x, y), point) || (alaskaPoint.point(x, y), point) || (hawaiiPoint.point(x, y), point);
  }
  albersUsa.invert = function(coordinates) {
    var k2 = lower48.scale(), t = lower48.translate(), x = (coordinates[0] - t[0]) / k2, y = (coordinates[1] - t[1]) / k2;
    return (y >= 0.12 && y < 0.234 && x >= -0.425 && x < -0.214 ? alaska : y >= 0.166 && y < 0.234 && x >= -0.214 && x < -0.115 ? hawaii : lower48).invert(coordinates);
  };
  albersUsa.stream = function(stream) {
    return cache && cacheStream === stream ? cache : cache = multiplex([lower48.stream(cacheStream = stream), alaska.stream(stream), hawaii.stream(stream)]);
  };
  albersUsa.precision = function(_) {
    if (!arguments.length)
      return lower48.precision();
    lower48.precision(_), alaska.precision(_), hawaii.precision(_);
    return reset();
  };
  albersUsa.scale = function(_) {
    if (!arguments.length)
      return lower48.scale();
    lower48.scale(_), alaska.scale(_ * 0.35), hawaii.scale(_);
    return albersUsa.translate(lower48.translate());
  };
  albersUsa.translate = function(_) {
    if (!arguments.length)
      return lower48.translate();
    var k2 = lower48.scale(), x = +_[0], y = +_[1];
    lower48Point = lower48.translate(_).clipExtent([[x - 0.455 * k2, y - 0.238 * k2], [x + 0.455 * k2, y + 0.238 * k2]]).stream(pointStream);
    alaskaPoint = alaska.translate([x - 0.307 * k2, y + 0.201 * k2]).clipExtent([[x - 0.425 * k2 + epsilon, y + 0.12 * k2 + epsilon], [x - 0.214 * k2 - epsilon, y + 0.234 * k2 - epsilon]]).stream(pointStream);
    hawaiiPoint = hawaii.translate([x - 0.205 * k2, y + 0.212 * k2]).clipExtent([[x - 0.214 * k2 + epsilon, y + 0.166 * k2 + epsilon], [x - 0.115 * k2 - epsilon, y + 0.234 * k2 - epsilon]]).stream(pointStream);
    return reset();
  };
  albersUsa.fitExtent = function(extent, object) {
    return fitExtent(albersUsa, extent, object);
  };
  albersUsa.fitSize = function(size, object) {
    return fitSize(albersUsa, size, object);
  };
  albersUsa.fitWidth = function(width, object) {
    return fitWidth(albersUsa, width, object);
  };
  albersUsa.fitHeight = function(height, object) {
    return fitHeight(albersUsa, height, object);
  };
  function reset() {
    cache = cacheStream = null;
    return albersUsa;
  }
  return albersUsa.scale(1070);
}

// node_modules/d3-geo/src/projection/azimuthal.js
function azimuthalRaw(scale) {
  return function(x, y) {
    var cx = cos(x), cy = cos(y), k2 = scale(cx * cy);
    if (k2 === Infinity)
      return [2, 0];
    return [
      k2 * cy * sin(x),
      k2 * sin(y)
    ];
  };
}
function azimuthalInvert(angle3) {
  return function(x, y) {
    var z = sqrt(x * x + y * y), c = angle3(z), sc = sin(c), cc = cos(c);
    return [
      atan2(x * sc, z * cc),
      asin(z && y * sc / z)
    ];
  };
}

// node_modules/d3-geo/src/projection/azimuthalEqualArea.js
var azimuthalEqualAreaRaw = azimuthalRaw(function(cxcy) {
  return sqrt(2 / (1 + cxcy));
});
azimuthalEqualAreaRaw.invert = azimuthalInvert(function(z) {
  return 2 * asin(z / 2);
});
function azimuthalEqualArea_default() {
  return projection(azimuthalEqualAreaRaw).scale(124.75).clipAngle(180 - 1e-3);
}

// node_modules/d3-geo/src/projection/azimuthalEquidistant.js
var azimuthalEquidistantRaw = azimuthalRaw(function(c) {
  return (c = acos(c)) && c / sin(c);
});
azimuthalEquidistantRaw.invert = azimuthalInvert(function(z) {
  return z;
});
function azimuthalEquidistant_default() {
  return projection(azimuthalEquidistantRaw).scale(79.4188).clipAngle(180 - 1e-3);
}

// node_modules/d3-geo/src/projection/mercator.js
function mercatorRaw(lambda, phi) {
  return [lambda, log(tan((halfPi + phi) / 2))];
}
mercatorRaw.invert = function(x, y) {
  return [x, 2 * atan(exp(y)) - halfPi];
};
function mercator_default() {
  return mercatorProjection(mercatorRaw).scale(961 / tau);
}
function mercatorProjection(project) {
  var m = projection(project), center = m.center, scale = m.scale, translate = m.translate, clipExtent = m.clipExtent, x03 = null, y03, x12, y12;
  m.scale = function(_) {
    return arguments.length ? (scale(_), reclip()) : scale();
  };
  m.translate = function(_) {
    return arguments.length ? (translate(_), reclip()) : translate();
  };
  m.center = function(_) {
    return arguments.length ? (center(_), reclip()) : center();
  };
  m.clipExtent = function(_) {
    return arguments.length ? (_ == null ? x03 = y03 = x12 = y12 = null : (x03 = +_[0][0], y03 = +_[0][1], x12 = +_[1][0], y12 = +_[1][1]), reclip()) : x03 == null ? null : [[x03, y03], [x12, y12]];
  };
  function reclip() {
    var k2 = pi * scale(), t = m(rotation_default(m.rotate()).invert([0, 0]));
    return clipExtent(x03 == null ? [[t[0] - k2, t[1] - k2], [t[0] + k2, t[1] + k2]] : project === mercatorRaw ? [[Math.max(t[0] - k2, x03), y03], [Math.min(t[0] + k2, x12), y12]] : [[x03, Math.max(t[1] - k2, y03)], [x12, Math.min(t[1] + k2, y12)]]);
  }
  return reclip();
}

// node_modules/d3-geo/src/projection/conicConformal.js
function tany(y) {
  return tan((halfPi + y) / 2);
}
function conicConformalRaw(y03, y12) {
  var cy0 = cos(y03), n = y03 === y12 ? sin(y03) : log(cy0 / cos(y12)) / log(tany(y12) / tany(y03)), f = cy0 * pow(tany(y03), n) / n;
  if (!n)
    return mercatorRaw;
  function project(x, y) {
    if (f > 0) {
      if (y < -halfPi + epsilon)
        y = -halfPi + epsilon;
    } else {
      if (y > halfPi - epsilon)
        y = halfPi - epsilon;
    }
    var r = f / pow(tany(y), n);
    return [r * sin(n * x), f - r * cos(n * x)];
  }
  project.invert = function(x, y) {
    var fy = f - y, r = sign(n) * sqrt(x * x + fy * fy), l = atan2(x, abs(fy)) * sign(fy);
    if (fy * n < 0)
      l -= pi * sign(x) * sign(fy);
    return [l / n, 2 * atan(pow(f / r, 1 / n)) - halfPi];
  };
  return project;
}
function conicConformal_default() {
  return conicProjection(conicConformalRaw).scale(109.5).parallels([30, 30]);
}

// node_modules/d3-geo/src/projection/equirectangular.js
function equirectangularRaw(lambda, phi) {
  return [lambda, phi];
}
equirectangularRaw.invert = equirectangularRaw;
function equirectangular_default() {
  return projection(equirectangularRaw).scale(152.63);
}

// node_modules/d3-geo/src/projection/conicEquidistant.js
function conicEquidistantRaw(y03, y12) {
  var cy0 = cos(y03), n = y03 === y12 ? sin(y03) : (cy0 - cos(y12)) / (y12 - y03), g = cy0 / n + y03;
  if (abs(n) < epsilon)
    return equirectangularRaw;
  function project(x, y) {
    var gy = g - y, nx = n * x;
    return [gy * sin(nx), g - gy * cos(nx)];
  }
  project.invert = function(x, y) {
    var gy = g - y, l = atan2(x, abs(gy)) * sign(gy);
    if (gy * n < 0)
      l -= pi * sign(x) * sign(gy);
    return [l / n, g - sign(n) * sqrt(x * x + gy * gy)];
  };
  return project;
}
function conicEquidistant_default() {
  return conicProjection(conicEquidistantRaw).scale(131.154).center([0, 13.9389]);
}

// node_modules/d3-geo/src/projection/equalEarth.js
var A1 = 1.340264;
var A2 = -0.081106;
var A3 = 893e-6;
var A4 = 3796e-6;
var M = sqrt(3) / 2;
var iterations = 12;
function equalEarthRaw(lambda, phi) {
  var l = asin(M * sin(phi)), l2 = l * l, l6 = l2 * l2 * l2;
  return [
    lambda * cos(l) / (M * (A1 + 3 * A2 * l2 + l6 * (7 * A3 + 9 * A4 * l2))),
    l * (A1 + A2 * l2 + l6 * (A3 + A4 * l2))
  ];
}
equalEarthRaw.invert = function(x, y) {
  var l = y, l2 = l * l, l6 = l2 * l2 * l2;
  for (var i = 0, delta, fy, fpy; i < iterations; ++i) {
    fy = l * (A1 + A2 * l2 + l6 * (A3 + A4 * l2)) - y;
    fpy = A1 + 3 * A2 * l2 + l6 * (7 * A3 + 9 * A4 * l2);
    l -= delta = fy / fpy, l2 = l * l, l6 = l2 * l2 * l2;
    if (abs(delta) < epsilon2)
      break;
  }
  return [
    M * x * (A1 + 3 * A2 * l2 + l6 * (7 * A3 + 9 * A4 * l2)) / cos(l),
    asin(sin(l) / M)
  ];
};
function equalEarth_default() {
  return projection(equalEarthRaw).scale(177.158);
}

// node_modules/d3-geo/src/projection/gnomonic.js
function gnomonicRaw(x, y) {
  var cy = cos(y), k2 = cos(x) * cy;
  return [cy * sin(x) / k2, sin(y) / k2];
}
gnomonicRaw.invert = azimuthalInvert(atan);
function gnomonic_default() {
  return projection(gnomonicRaw).scale(144.049).clipAngle(60);
}

// node_modules/d3-geo/src/projection/naturalEarth1.js
function naturalEarth1Raw(lambda, phi) {
  var phi2 = phi * phi, phi4 = phi2 * phi2;
  return [
    lambda * (0.8707 - 0.131979 * phi2 + phi4 * (-0.013791 + phi4 * (3971e-6 * phi2 - 1529e-6 * phi4))),
    phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 5916e-6 * phi4)))
  ];
}
naturalEarth1Raw.invert = function(x, y) {
  var phi = y, i = 25, delta;
  do {
    var phi2 = phi * phi, phi4 = phi2 * phi2;
    phi -= delta = (phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 5916e-6 * phi4))) - y) / (1.007226 + phi2 * (0.015085 * 3 + phi4 * (-0.044475 * 7 + 0.028874 * 9 * phi2 - 5916e-6 * 11 * phi4)));
  } while (abs(delta) > epsilon && --i > 0);
  return [
    x / (0.8707 + (phi2 = phi * phi) * (-0.131979 + phi2 * (-0.013791 + phi2 * phi2 * phi2 * (3971e-6 - 1529e-6 * phi2)))),
    phi
  ];
};
function naturalEarth1_default() {
  return projection(naturalEarth1Raw).scale(175.295);
}

// node_modules/d3-geo/src/projection/orthographic.js
function orthographicRaw(x, y) {
  return [cos(y) * sin(x), sin(y)];
}
orthographicRaw.invert = azimuthalInvert(asin);
function orthographic_default() {
  return projection(orthographicRaw).scale(249.5).clipAngle(90 + epsilon);
}

// node_modules/d3-geo/src/projection/stereographic.js
function stereographicRaw(x, y) {
  var cy = cos(y), k2 = 1 + cos(x) * cy;
  return [cy * sin(x) / k2, sin(y) / k2];
}
stereographicRaw.invert = azimuthalInvert(function(z) {
  return 2 * atan(z);
});
function stereographic_default() {
  return projection(stereographicRaw).scale(250).clipAngle(142);
}

// node_modules/d3-geo/src/projection/transverseMercator.js
function transverseMercatorRaw(lambda, phi) {
  return [log(tan((halfPi + phi) / 2)), -lambda];
}
transverseMercatorRaw.invert = function(x, y) {
  return [-y, 2 * atan(exp(x)) - halfPi];
};
function transverseMercator_default() {
  var m = mercatorProjection(transverseMercatorRaw), center = m.center, rotate = m.rotate;
  m.center = function(_) {
    return arguments.length ? center([-_[1], _[0]]) : (_ = center(), [_[1], -_[0]]);
  };
  m.rotate = function(_) {
    return arguments.length ? rotate([_[0], _[1], _.length > 2 ? _[2] + 90 : 90]) : (_ = rotate(), [_[0], _[1], _[2] - 90]);
  };
  return rotate([0, 0, 90]).scale(159.155);
}

// node_modules/d3-geo-projection/src/math.js
var abs2 = Math.abs;
var atan3 = Math.atan;
var atan22 = Math.atan2;
var cos2 = Math.cos;
var exp2 = Math.exp;
var floor = Math.floor;
var log2 = Math.log;
var max = Math.max;
var min = Math.min;
var pow2 = Math.pow;
var round = Math.round;
var sign2 = Math.sign || function(x) {
  return x > 0 ? 1 : x < 0 ? -1 : 0;
};
var sin2 = Math.sin;
var tan2 = Math.tan;
var epsilon3 = 1e-6;
var epsilon22 = 1e-12;
var pi2 = Math.PI;
var halfPi2 = pi2 / 2;
var quarterPi2 = pi2 / 4;
var sqrt1_2 = Math.SQRT1_2;
var sqrt2 = sqrt3(2);
var sqrtPi = sqrt3(pi2);
var tau2 = pi2 * 2;
var degrees2 = 180 / pi2;
var radians2 = pi2 / 180;
function sinci(x) {
  return x ? x / Math.sin(x) : 1;
}
function asin2(x) {
  return x > 1 ? halfPi2 : x < -1 ? -halfPi2 : Math.asin(x);
}
function acos2(x) {
  return x > 1 ? 0 : x < -1 ? pi2 : Math.acos(x);
}
function sqrt3(x) {
  return x > 0 ? Math.sqrt(x) : 0;
}
function tanh(x) {
  x = exp2(2 * x);
  return (x - 1) / (x + 1);
}
function sinh(x) {
  return (exp2(x) - exp2(-x)) / 2;
}
function cosh(x) {
  return (exp2(x) + exp2(-x)) / 2;
}
function arsinh(x) {
  return log2(x + sqrt3(x * x + 1));
}
function arcosh(x) {
  return log2(x + sqrt3(x * x - 1));
}

// node_modules/d3-geo-projection/src/airy.js
function airyRaw(beta) {
  var tanBeta_2 = tan2(beta / 2), b = 2 * log2(cos2(beta / 2)) / (tanBeta_2 * tanBeta_2);
  function forward(x, y) {
    var cosx = cos2(x), cosy = cos2(y), siny = sin2(y), cosz = cosy * cosx, k2 = -((1 - cosz ? log2((1 + cosz) / 2) / (1 - cosz) : -0.5) + b / (1 + cosz));
    return [k2 * cosy * sin2(x), k2 * siny];
  }
  forward.invert = function(x, y) {
    var r = sqrt3(x * x + y * y), z = -beta / 2, i = 50, delta;
    if (!r)
      return [0, 0];
    do {
      var z_2 = z / 2, cosz_2 = cos2(z_2), sinz_2 = sin2(z_2), tanz_2 = sinz_2 / cosz_2, lnsecz_2 = -log2(abs2(cosz_2));
      z -= delta = (2 / tanz_2 * lnsecz_2 - b * tanz_2 - r) / (-lnsecz_2 / (sinz_2 * sinz_2) + 1 - b / (2 * cosz_2 * cosz_2)) * (cosz_2 < 0 ? 0.7 : 1);
    } while (abs2(delta) > epsilon3 && --i > 0);
    var sinz = sin2(z);
    return [atan22(x * sinz, r * cos2(z)), asin2(y * sinz / r)];
  };
  return forward;
}
function airy_default() {
  var beta = halfPi2, m = projectionMutator(airyRaw), p = m(beta);
  p.radius = function(_) {
    return arguments.length ? m(beta = _ * radians2) : beta * degrees2;
  };
  return p.scale(179.976).clipAngle(147);
}

// node_modules/d3-geo-projection/src/aitoff.js
function aitoffRaw(x, y) {
  var cosy = cos2(y), sincia = sinci(acos2(cosy * cos2(x /= 2)));
  return [2 * cosy * sin2(x) * sincia, sin2(y) * sincia];
}
aitoffRaw.invert = function(x, y) {
  if (x * x + 4 * y * y > pi2 * pi2 + epsilon3)
    return;
  var x12 = x, y12 = y, i = 25;
  do {
    var sinx = sin2(x12), sinx_2 = sin2(x12 / 2), cosx_2 = cos2(x12 / 2), siny = sin2(y12), cosy = cos2(y12), sin_2y = sin2(2 * y12), sin2y = siny * siny, cos2y = cosy * cosy, sin2x_2 = sinx_2 * sinx_2, c = 1 - cos2y * cosx_2 * cosx_2, e = c ? acos2(cosy * cosx_2) * sqrt3(f = 1 / c) : f = 0, f, fx = 2 * e * cosy * sinx_2 - x, fy = e * siny - y, dxdx = f * (cos2y * sin2x_2 + e * cosy * cosx_2 * sin2y), dxdy = f * (0.5 * sinx * sin_2y - e * 2 * siny * sinx_2), dydx = f * 0.25 * (sin_2y * sinx_2 - e * siny * cos2y * sinx), dydy = f * (sin2y * cosx_2 + e * sin2x_2 * cosy), z = dxdy * dydx - dydy * dxdx;
    if (!z)
      break;
    var dx = (fy * dxdy - fx * dydy) / z, dy = (fx * dydx - fy * dxdx) / z;
    x12 -= dx, y12 -= dy;
  } while ((abs2(dx) > epsilon3 || abs2(dy) > epsilon3) && --i > 0);
  return [x12, y12];
};
function aitoff_default() {
  return projection(aitoffRaw).scale(152.63);
}

// node_modules/d3-geo-projection/src/armadillo.js
function armadilloRaw(phi03) {
  var sinPhi02 = sin2(phi03), cosPhi02 = cos2(phi03), sPhi0 = phi03 >= 0 ? 1 : -1, tanPhi0 = tan2(sPhi0 * phi03), k2 = (1 + sinPhi02 - cosPhi02) / 2;
  function forward(lambda, phi) {
    var cosPhi = cos2(phi), cosLambda = cos2(lambda /= 2);
    return [
      (1 + cosPhi) * sin2(lambda),
      (sPhi0 * phi > -atan22(cosLambda, tanPhi0) - 1e-3 ? 0 : -sPhi0 * 10) + k2 + sin2(phi) * cosPhi02 - (1 + cosPhi) * sinPhi02 * cosLambda
      // TODO D3 core should allow null or [NaN, NaN] to be returned.
    ];
  }
  forward.invert = function(x, y) {
    var lambda = 0, phi = 0, i = 50;
    do {
      var cosLambda = cos2(lambda), sinLambda = sin2(lambda), cosPhi = cos2(phi), sinPhi = sin2(phi), A5 = 1 + cosPhi, fx = A5 * sinLambda - x, fy = k2 + sinPhi * cosPhi02 - A5 * sinPhi02 * cosLambda - y, dxdLambda = A5 * cosLambda / 2, dxdPhi = -sinLambda * sinPhi, dydLambda = sinPhi02 * A5 * sinLambda / 2, dydPhi = cosPhi02 * cosPhi + sinPhi02 * cosLambda * sinPhi, denominator = dxdPhi * dydLambda - dydPhi * dxdLambda, dLambda = (fy * dxdPhi - fx * dydPhi) / denominator / 2, dPhi = (fx * dydLambda - fy * dxdLambda) / denominator;
      if (abs2(dPhi) > 2)
        dPhi /= 2;
      lambda -= dLambda, phi -= dPhi;
    } while ((abs2(dLambda) > epsilon3 || abs2(dPhi) > epsilon3) && --i > 0);
    return sPhi0 * phi > -atan22(cos2(lambda), tanPhi0) - 1e-3 ? [lambda * 2, phi] : null;
  };
  return forward;
}
function armadillo_default() {
  var phi03 = 20 * radians2, sPhi0 = phi03 >= 0 ? 1 : -1, tanPhi0 = tan2(sPhi0 * phi03), m = projectionMutator(armadilloRaw), p = m(phi03), stream_ = p.stream;
  p.parallel = function(_) {
    if (!arguments.length)
      return phi03 * degrees2;
    tanPhi0 = tan2((sPhi0 = (phi03 = _ * radians2) >= 0 ? 1 : -1) * phi03);
    return m(phi03);
  };
  p.stream = function(stream) {
    var rotate = p.rotate(), rotateStream = stream_(stream), sphereStream = (p.rotate([0, 0]), stream_(stream)), precision = p.precision();
    p.rotate(rotate);
    rotateStream.sphere = function() {
      sphereStream.polygonStart(), sphereStream.lineStart();
      for (var lambda = sPhi0 * -180; sPhi0 * lambda < 180; lambda += sPhi0 * 90)
        sphereStream.point(lambda, sPhi0 * 90);
      if (phi03)
        while (sPhi0 * (lambda -= 3 * sPhi0 * precision) >= -180) {
          sphereStream.point(lambda, sPhi0 * -atan22(cos2(lambda * radians2 / 2), tanPhi0) * degrees2);
        }
      sphereStream.lineEnd(), sphereStream.polygonEnd();
    };
    return rotateStream;
  };
  return p.scale(218.695).center([0, 28.0974]);
}

// node_modules/d3-geo-projection/src/august.js
function augustRaw(lambda, phi) {
  var tanPhi = tan2(phi / 2), k2 = sqrt3(1 - tanPhi * tanPhi), c = 1 + k2 * cos2(lambda /= 2), x = sin2(lambda) * k2 / c, y = tanPhi / c, x2 = x * x, y2 = y * y;
  return [
    4 / 3 * x * (3 + x2 - 3 * y2),
    4 / 3 * y * (3 + 3 * x2 - y2)
  ];
}
augustRaw.invert = function(x, y) {
  x *= 3 / 8, y *= 3 / 8;
  if (!x && abs2(y) > 1)
    return null;
  var x2 = x * x, y2 = y * y, s = 1 + x2 + y2, sin3Eta = sqrt3((s - sqrt3(s * s - 4 * y * y)) / 2), eta = asin2(sin3Eta) / 3, xi = sin3Eta ? arcosh(abs2(y / sin3Eta)) / 3 : arsinh(abs2(x)) / 3, cosEta = cos2(eta), coshXi = cosh(xi), d = coshXi * coshXi - cosEta * cosEta;
  return [
    sign2(x) * 2 * atan22(sinh(xi) * cosEta, 0.25 - d),
    sign2(y) * 2 * atan22(coshXi * sin2(eta), 0.25 + d)
  ];
};
function august_default() {
  return projection(augustRaw).scale(66.1603);
}

// node_modules/d3-geo-projection/src/baker.js
var sqrt8 = sqrt3(8);
var phi02 = log2(1 + sqrt2);
function bakerRaw(lambda, phi) {
  var phi03 = abs2(phi);
  return phi03 < quarterPi2 ? [lambda, log2(tan2(quarterPi2 + phi / 2))] : [lambda * cos2(phi03) * (2 * sqrt2 - 1 / sin2(phi03)), sign2(phi) * (2 * sqrt2 * (phi03 - quarterPi2) - log2(tan2(phi03 / 2)))];
}
bakerRaw.invert = function(x, y) {
  if ((y03 = abs2(y)) < phi02)
    return [x, 2 * atan3(exp2(y)) - halfPi2];
  var phi = quarterPi2, i = 25, delta, y03;
  do {
    var cosPhi_2 = cos2(phi / 2), tanPhi_2 = tan2(phi / 2);
    phi -= delta = (sqrt8 * (phi - quarterPi2) - log2(tanPhi_2) - y03) / (sqrt8 - cosPhi_2 * cosPhi_2 / (2 * tanPhi_2));
  } while (abs2(delta) > epsilon22 && --i > 0);
  return [x / (cos2(phi) * (sqrt8 - 1 / sin2(phi))), sign2(y) * phi];
};
function baker_default() {
  return projection(bakerRaw).scale(112.314);
}

// node_modules/d3-geo-projection/src/berghaus.js
function berghausRaw(lobes7) {
  var k2 = 2 * pi2 / lobes7;
  function forward(lambda, phi) {
    var p = azimuthalEquidistantRaw(lambda, phi);
    if (abs2(lambda) > halfPi2) {
      var theta = atan22(p[1], p[0]), r = sqrt3(p[0] * p[0] + p[1] * p[1]), theta0 = k2 * round((theta - halfPi2) / k2) + halfPi2, alpha = atan22(sin2(theta -= theta0), 2 - cos2(theta));
      theta = theta0 + asin2(pi2 / r * sin2(alpha)) - alpha;
      p[0] = r * cos2(theta);
      p[1] = r * sin2(theta);
    }
    return p;
  }
  forward.invert = function(x, y) {
    var r = sqrt3(x * x + y * y);
    if (r > halfPi2) {
      var theta = atan22(y, x), theta0 = k2 * round((theta - halfPi2) / k2) + halfPi2, s = theta > theta0 ? -1 : 1, A5 = r * cos2(theta0 - theta), cotAlpha = 1 / tan2(s * acos2((A5 - pi2) / sqrt3(pi2 * (pi2 - 2 * A5) + r * r)));
      theta = theta0 + 2 * atan3((cotAlpha + s * sqrt3(cotAlpha * cotAlpha - 3)) / 3);
      x = r * cos2(theta), y = r * sin2(theta);
    }
    return azimuthalEquidistantRaw.invert(x, y);
  };
  return forward;
}
function berghaus_default() {
  var lobes7 = 5, m = projectionMutator(berghausRaw), p = m(lobes7), projectionStream = p.stream, epsilon4 = 0.01, cr = -cos2(epsilon4 * radians2), sr = sin2(epsilon4 * radians2);
  p.lobes = function(_) {
    return arguments.length ? m(lobes7 = +_) : lobes7;
  };
  p.stream = function(stream) {
    var rotate = p.rotate(), rotateStream = projectionStream(stream), sphereStream = (p.rotate([0, 0]), projectionStream(stream));
    p.rotate(rotate);
    rotateStream.sphere = function() {
      sphereStream.polygonStart(), sphereStream.lineStart();
      for (var i = 0, delta = 360 / lobes7, delta0 = 2 * pi2 / lobes7, phi = 90 - 180 / lobes7, phi03 = halfPi2; i < lobes7; ++i, phi -= delta, phi03 -= delta0) {
        sphereStream.point(atan22(sr * cos2(phi03), cr) * degrees2, asin2(sr * sin2(phi03)) * degrees2);
        if (phi < -90) {
          sphereStream.point(-90, -180 - phi - epsilon4);
          sphereStream.point(-90, -180 - phi + epsilon4);
        } else {
          sphereStream.point(90, phi + epsilon4);
          sphereStream.point(90, phi - epsilon4);
        }
      }
      sphereStream.lineEnd(), sphereStream.polygonEnd();
    };
    return rotateStream;
  };
  return p.scale(87.8076).center([0, 17.1875]).clipAngle(180 - 1e-3);
}

// node_modules/d3-geo-projection/src/hammer.js
function hammerRaw(A5, B2) {
  if (arguments.length < 2)
    B2 = A5;
  if (B2 === 1)
    return azimuthalEqualAreaRaw;
  if (B2 === Infinity)
    return hammerQuarticAuthalicRaw;
  function forward(lambda, phi) {
    var coordinates = azimuthalEqualAreaRaw(lambda / B2, phi);
    coordinates[0] *= A5;
    return coordinates;
  }
  forward.invert = function(x, y) {
    var coordinates = azimuthalEqualAreaRaw.invert(x / A5, y);
    coordinates[0] *= B2;
    return coordinates;
  };
  return forward;
}
function hammerQuarticAuthalicRaw(lambda, phi) {
  return [
    lambda * cos2(phi) / cos2(phi /= 2),
    2 * sin2(phi)
  ];
}
hammerQuarticAuthalicRaw.invert = function(x, y) {
  var phi = 2 * asin2(y / 2);
  return [
    x * cos2(phi / 2) / cos2(phi),
    phi
  ];
};
function hammer_default() {
  var B2 = 2, m = projectionMutator(hammerRaw), p = m(B2);
  p.coefficient = function(_) {
    if (!arguments.length)
      return B2;
    return m(B2 = +_);
  };
  return p.scale(169.529);
}

// node_modules/d3-geo-projection/src/newton.js
function solve2d(f, MAX_ITERATIONS, eps) {
  if (MAX_ITERATIONS === void 0)
    MAX_ITERATIONS = 40;
  if (eps === void 0)
    eps = epsilon22;
  return function(x, y, a, b) {
    var err2, da, db;
    a = a === void 0 ? 0 : +a;
    b = b === void 0 ? 0 : +b;
    for (var i = 0; i < MAX_ITERATIONS; i++) {
      var p = f(a, b), tx = p[0] - x, ty = p[1] - y;
      if (abs2(tx) < eps && abs2(ty) < eps)
        break;
      var h = tx * tx + ty * ty;
      if (h > err2) {
        a -= da /= 2;
        b -= db /= 2;
        continue;
      }
      err2 = h;
      var ea = (a > 0 ? -1 : 1) * eps, eb = (b > 0 ? -1 : 1) * eps, pa = f(a + ea, b), pb = f(a, b + eb), dxa = (pa[0] - p[0]) / ea, dya = (pa[1] - p[1]) / ea, dxb = (pb[0] - p[0]) / eb, dyb = (pb[1] - p[1]) / eb, D = dyb * dxa - dya * dxb, l = (abs2(D) < 0.5 ? 0.5 : 1) / D;
      da = (ty * dxb - tx * dyb) * l;
      db = (tx * dya - ty * dxa) * l;
      a += da;
      b += db;
      if (abs2(da) < eps && abs2(db) < eps)
        break;
    }
    return [a, b];
  };
}

// node_modules/d3-geo-projection/src/mollweide.js
function mollweideBromleyTheta(cp, phi) {
  var cpsinPhi = cp * sin2(phi), i = 30, delta;
  do
    phi -= delta = (phi + sin2(phi) - cpsinPhi) / (1 + cos2(phi));
  while (abs2(delta) > epsilon3 && --i > 0);
  return phi / 2;
}
function mollweideBromleyRaw(cx, cy, cp) {
  function forward(lambda, phi) {
    return [cx * lambda * cos2(phi = mollweideBromleyTheta(cp, phi)), cy * sin2(phi)];
  }
  forward.invert = function(x, y) {
    return y = asin2(y / cy), [x / (cx * cos2(y)), asin2((2 * y + sin2(2 * y)) / cp)];
  };
  return forward;
}
var mollweideRaw = mollweideBromleyRaw(sqrt2 / halfPi2, sqrt2, pi2);
function mollweide_default() {
  return projection(mollweideRaw).scale(169.529);
}

// node_modules/d3-geo-projection/src/boggs.js
var k = 2.00276;
var w = 1.11072;
function boggsRaw(lambda, phi) {
  var theta = mollweideBromleyTheta(pi2, phi);
  return [k * lambda / (1 / cos2(phi) + w / cos2(theta)), (phi + sqrt2 * sin2(theta)) / k];
}
boggsRaw.invert = function(x, y) {
  var ky = k * y, theta = y < 0 ? -quarterPi2 : quarterPi2, i = 25, delta, phi;
  do {
    phi = ky - sqrt2 * sin2(theta);
    theta -= delta = (sin2(2 * theta) + 2 * theta - pi2 * sin2(phi)) / (2 * cos2(2 * theta) + 2 + pi2 * cos2(phi) * sqrt2 * cos2(theta));
  } while (abs2(delta) > epsilon3 && --i > 0);
  phi = ky - sqrt2 * sin2(theta);
  return [x * (1 / cos2(phi) + w / cos2(theta)) / k, phi];
};
function boggs_default() {
  return projection(boggsRaw).scale(160.857);
}

// node_modules/d3-geo-projection/src/parallel1.js
function parallel1_default(projectAt) {
  var phi03 = 0, m = projectionMutator(projectAt), p = m(phi03);
  p.parallel = function(_) {
    return arguments.length ? m(phi03 = _ * radians2) : phi03 * degrees2;
  };
  return p;
}

// node_modules/d3-geo-projection/src/sinusoidal.js
function sinusoidalRaw(lambda, phi) {
  return [lambda * cos2(phi), phi];
}
sinusoidalRaw.invert = function(x, y) {
  return [x / cos2(y), y];
};
function sinusoidal_default() {
  return projection(sinusoidalRaw).scale(152.63);
}

// node_modules/d3-geo-projection/src/bonne.js
function bonneRaw(phi03) {
  if (!phi03)
    return sinusoidalRaw;
  var cotPhi0 = 1 / tan2(phi03);
  function forward(lambda, phi) {
    var rho = cotPhi0 + phi03 - phi, e = rho ? lambda * cos2(phi) / rho : rho;
    return [rho * sin2(e), cotPhi0 - rho * cos2(e)];
  }
  forward.invert = function(x, y) {
    var rho = sqrt3(x * x + (y = cotPhi0 - y) * y), phi = cotPhi0 + phi03 - rho;
    return [rho / cos2(phi) * atan22(x, y), phi];
  };
  return forward;
}
function bonne_default() {
  return parallel1_default(bonneRaw).scale(123.082).center([0, 26.1441]).parallel(45);
}

// node_modules/d3-geo-projection/src/bottomley.js
function bottomleyRaw(sinPsi) {
  function forward(lambda, phi) {
    var rho = halfPi2 - phi, eta = rho ? lambda * sinPsi * sin2(rho) / rho : rho;
    return [rho * sin2(eta) / sinPsi, halfPi2 - rho * cos2(eta)];
  }
  forward.invert = function(x, y) {
    var x12 = x * sinPsi, y12 = halfPi2 - y, rho = sqrt3(x12 * x12 + y12 * y12), eta = atan22(x12, y12);
    return [(rho ? rho / sin2(rho) : 1) * eta / sinPsi, halfPi2 - rho];
  };
  return forward;
}
function bottomley_default() {
  var sinPsi = 0.5, m = projectionMutator(bottomleyRaw), p = m(sinPsi);
  p.fraction = function(_) {
    return arguments.length ? m(sinPsi = +_) : sinPsi;
  };
  return p.scale(158.837);
}

// node_modules/d3-geo-projection/src/bromley.js
var bromleyRaw = mollweideBromleyRaw(1, 4 / pi2, pi2);
function bromley_default() {
  return projection(bromleyRaw).scale(152.63);
}

// node_modules/d3-geo-projection/src/collignon.js
function collignonRaw(lambda, phi) {
  var alpha = sqrt3(1 - sin2(phi));
  return [2 / sqrtPi * lambda * alpha, sqrtPi * (1 - alpha)];
}
collignonRaw.invert = function(x, y) {
  var lambda = (lambda = y / sqrtPi - 1) * lambda;
  return [lambda > 0 ? x * sqrt3(pi2 / lambda) / 2 : 0, asin2(1 - lambda)];
};
function collignon_default() {
  return projection(collignonRaw).scale(95.6464).center([0, 30]);
}

// node_modules/d3-geo-projection/src/craig.js
function craigRaw(phi03) {
  var tanPhi0 = tan2(phi03);
  function forward(lambda, phi) {
    return [lambda, (lambda ? lambda / sin2(lambda) : 1) * (sin2(phi) * cos2(lambda) - tanPhi0 * cos2(phi))];
  }
  forward.invert = tanPhi0 ? function(x, y) {
    if (x)
      y *= sin2(x) / x;
    var cosLambda = cos2(x);
    return [x, 2 * atan22(sqrt3(cosLambda * cosLambda + tanPhi0 * tanPhi0 - y * y) - cosLambda, tanPhi0 - y)];
  } : function(x, y) {
    return [x, asin2(x ? y * tan2(x) / x : y)];
  };
  return forward;
}
function craig_default() {
  return parallel1_default(craigRaw).scale(249.828).clipAngle(90);
}

// node_modules/d3-geo-projection/src/craster.js
var sqrt32 = sqrt3(3);
function crasterRaw(lambda, phi) {
  return [sqrt32 * lambda * (2 * cos2(2 * phi / 3) - 1) / sqrtPi, sqrt32 * sqrtPi * sin2(phi / 3)];
}
crasterRaw.invert = function(x, y) {
  var phi = 3 * asin2(y / (sqrt32 * sqrtPi));
  return [sqrtPi * x / (sqrt32 * (2 * cos2(2 * phi / 3) - 1)), phi];
};
function craster_default() {
  return projection(crasterRaw).scale(156.19);
}

// node_modules/d3-geo-projection/src/cylindricalEqualArea.js
function cylindricalEqualAreaRaw2(phi03) {
  var cosPhi02 = cos2(phi03);
  function forward(lambda, phi) {
    return [lambda * cosPhi02, sin2(phi) / cosPhi02];
  }
  forward.invert = function(x, y) {
    return [x / cosPhi02, asin2(y * cosPhi02)];
  };
  return forward;
}
function cylindricalEqualArea_default() {
  return parallel1_default(cylindricalEqualAreaRaw2).parallel(38.58).scale(195.044);
}

// node_modules/d3-geo-projection/src/cylindricalStereographic.js
function cylindricalStereographicRaw(phi03) {
  var cosPhi02 = cos2(phi03);
  function forward(lambda, phi) {
    return [lambda * cosPhi02, (1 + cosPhi02) * tan2(phi / 2)];
  }
  forward.invert = function(x, y) {
    return [x / cosPhi02, atan3(y / (1 + cosPhi02)) * 2];
  };
  return forward;
}
function cylindricalStereographic_default() {
  return parallel1_default(cylindricalStereographicRaw).scale(124.75);
}

// node_modules/d3-geo-projection/src/eckert1.js
function eckert1Raw(lambda, phi) {
  var alpha = sqrt3(8 / (3 * pi2));
  return [
    alpha * lambda * (1 - abs2(phi) / pi2),
    alpha * phi
  ];
}
eckert1Raw.invert = function(x, y) {
  var alpha = sqrt3(8 / (3 * pi2)), phi = y / alpha;
  return [
    x / (alpha * (1 - abs2(phi) / pi2)),
    phi
  ];
};
function eckert1_default() {
  return projection(eckert1Raw).scale(165.664);
}

// node_modules/d3-geo-projection/src/eckert2.js
function eckert2Raw(lambda, phi) {
  var alpha = sqrt3(4 - 3 * sin2(abs2(phi)));
  return [
    2 / sqrt3(6 * pi2) * lambda * alpha,
    sign2(phi) * sqrt3(2 * pi2 / 3) * (2 - alpha)
  ];
}
eckert2Raw.invert = function(x, y) {
  var alpha = 2 - abs2(y) / sqrt3(2 * pi2 / 3);
  return [
    x * sqrt3(6 * pi2) / (2 * alpha),
    sign2(y) * asin2((4 - alpha * alpha) / 3)
  ];
};
function eckert2_default() {
  return projection(eckert2Raw).scale(165.664);
}

// node_modules/d3-geo-projection/src/eckert3.js
function eckert3Raw(lambda, phi) {
  var k2 = sqrt3(pi2 * (4 + pi2));
  return [
    2 / k2 * lambda * (1 + sqrt3(1 - 4 * phi * phi / (pi2 * pi2))),
    4 / k2 * phi
  ];
}
eckert3Raw.invert = function(x, y) {
  var k2 = sqrt3(pi2 * (4 + pi2)) / 2;
  return [
    x * k2 / (1 + sqrt3(1 - y * y * (4 + pi2) / (4 * pi2))),
    y * k2 / 2
  ];
};
function eckert3_default() {
  return projection(eckert3Raw).scale(180.739);
}

// node_modules/d3-geo-projection/src/eckert4.js
function eckert4Raw(lambda, phi) {
  var k2 = (2 + halfPi2) * sin2(phi);
  phi /= 2;
  for (var i = 0, delta = Infinity; i < 10 && abs2(delta) > epsilon3; i++) {
    var cosPhi = cos2(phi);
    phi -= delta = (phi + sin2(phi) * (cosPhi + 2) - k2) / (2 * cosPhi * (1 + cosPhi));
  }
  return [
    2 / sqrt3(pi2 * (4 + pi2)) * lambda * (1 + cos2(phi)),
    2 * sqrt3(pi2 / (4 + pi2)) * sin2(phi)
  ];
}
eckert4Raw.invert = function(x, y) {
  var A5 = y * sqrt3((4 + pi2) / pi2) / 2, k2 = asin2(A5), c = cos2(k2);
  return [
    x / (2 / sqrt3(pi2 * (4 + pi2)) * (1 + c)),
    asin2((k2 + A5 * (c + 2)) / (2 + halfPi2))
  ];
};
function eckert4_default() {
  return projection(eckert4Raw).scale(180.739);
}

// node_modules/d3-geo-projection/src/eckert5.js
function eckert5Raw(lambda, phi) {
  return [
    lambda * (1 + cos2(phi)) / sqrt3(2 + pi2),
    2 * phi / sqrt3(2 + pi2)
  ];
}
eckert5Raw.invert = function(x, y) {
  var k2 = sqrt3(2 + pi2), phi = y * k2 / 2;
  return [
    k2 * x / (1 + cos2(phi)),
    phi
  ];
};
function eckert5_default() {
  return projection(eckert5Raw).scale(173.044);
}

// node_modules/d3-geo-projection/src/eckert6.js
function eckert6Raw(lambda, phi) {
  var k2 = (1 + halfPi2) * sin2(phi);
  for (var i = 0, delta = Infinity; i < 10 && abs2(delta) > epsilon3; i++) {
    phi -= delta = (phi + sin2(phi) - k2) / (1 + cos2(phi));
  }
  k2 = sqrt3(2 + pi2);
  return [
    lambda * (1 + cos2(phi)) / k2,
    2 * phi / k2
  ];
}
eckert6Raw.invert = function(x, y) {
  var j = 1 + halfPi2, k2 = sqrt3(j / 2);
  return [
    x * 2 * k2 / (1 + cos2(y *= k2)),
    asin2((y + sin2(y)) / j)
  ];
};
function eckert6_default() {
  return projection(eckert6Raw).scale(173.044);
}

// node_modules/d3-geo-projection/src/eisenlohr.js
var eisenlohrK = 3 + 2 * sqrt2;
function eisenlohrRaw(lambda, phi) {
  var s0 = sin2(lambda /= 2), c0 = cos2(lambda), k2 = sqrt3(cos2(phi)), c1 = cos2(phi /= 2), t = sin2(phi) / (c1 + sqrt2 * c0 * k2), c = sqrt3(2 / (1 + t * t)), v = sqrt3((sqrt2 * c1 + (c0 + s0) * k2) / (sqrt2 * c1 + (c0 - s0) * k2));
  return [
    eisenlohrK * (c * (v - 1 / v) - 2 * log2(v)),
    eisenlohrK * (c * t * (v + 1 / v) - 2 * atan3(t))
  ];
}
eisenlohrRaw.invert = function(x, y) {
  if (!(p = augustRaw.invert(x / 1.2, y * 1.065)))
    return null;
  var lambda = p[0], phi = p[1], i = 20, p;
  x /= eisenlohrK, y /= eisenlohrK;
  do {
    var _0 = lambda / 2, _1 = phi / 2, s0 = sin2(_0), c0 = cos2(_0), s1 = sin2(_1), c1 = cos2(_1), cos1 = cos2(phi), k2 = sqrt3(cos1), t = s1 / (c1 + sqrt2 * c0 * k2), t2 = t * t, c = sqrt3(2 / (1 + t2)), v0 = sqrt2 * c1 + (c0 + s0) * k2, v1 = sqrt2 * c1 + (c0 - s0) * k2, v2 = v0 / v1, v = sqrt3(v2), vm1v = v - 1 / v, vp1v = v + 1 / v, fx = c * vm1v - 2 * log2(v) - x, fy = c * t * vp1v - 2 * atan3(t) - y, deltatDeltaLambda = s1 && sqrt1_2 * k2 * s0 * t2 / s1, deltatDeltaPhi = (sqrt2 * c0 * c1 + k2) / (2 * (c1 + sqrt2 * c0 * k2) * (c1 + sqrt2 * c0 * k2) * k2), deltacDeltat = -0.5 * t * c * c * c, deltacDeltaLambda = deltacDeltat * deltatDeltaLambda, deltacDeltaPhi = deltacDeltat * deltatDeltaPhi, A5 = (A5 = 2 * c1 + sqrt2 * k2 * (c0 - s0)) * A5 * v, deltavDeltaLambda = (sqrt2 * c0 * c1 * k2 + cos1) / A5, deltavDeltaPhi = -(sqrt2 * s0 * s1) / (k2 * A5), deltaxDeltaLambda = vm1v * deltacDeltaLambda - 2 * deltavDeltaLambda / v + c * (deltavDeltaLambda + deltavDeltaLambda / v2), deltaxDeltaPhi = vm1v * deltacDeltaPhi - 2 * deltavDeltaPhi / v + c * (deltavDeltaPhi + deltavDeltaPhi / v2), deltayDeltaLambda = t * vp1v * deltacDeltaLambda - 2 * deltatDeltaLambda / (1 + t2) + c * vp1v * deltatDeltaLambda + c * t * (deltavDeltaLambda - deltavDeltaLambda / v2), deltayDeltaPhi = t * vp1v * deltacDeltaPhi - 2 * deltatDeltaPhi / (1 + t2) + c * vp1v * deltatDeltaPhi + c * t * (deltavDeltaPhi - deltavDeltaPhi / v2), denominator = deltaxDeltaPhi * deltayDeltaLambda - deltayDeltaPhi * deltaxDeltaLambda;
    if (!denominator)
      break;
    var deltaLambda = (fy * deltaxDeltaPhi - fx * deltayDeltaPhi) / denominator, deltaPhi = (fx * deltayDeltaLambda - fy * deltaxDeltaLambda) / denominator;
    lambda -= deltaLambda;
    phi = max(-halfPi2, min(halfPi2, phi - deltaPhi));
  } while ((abs2(deltaLambda) > epsilon3 || abs2(deltaPhi) > epsilon3) && --i > 0);
  return abs2(abs2(phi) - halfPi2) < epsilon3 ? [0, phi] : i && [lambda, phi];
};
function eisenlohr_default() {
  return projection(eisenlohrRaw).scale(62.5271);
}

// node_modules/d3-geo-projection/src/fahey.js
var faheyK = cos2(35 * radians2);
function faheyRaw(lambda, phi) {
  var t = tan2(phi / 2);
  return [lambda * faheyK * sqrt3(1 - t * t), (1 + faheyK) * t];
}
faheyRaw.invert = function(x, y) {
  var t = y / (1 + faheyK);
  return [x && x / (faheyK * sqrt3(1 - t * t)), 2 * atan3(t)];
};
function fahey_default() {
  return projection(faheyRaw).scale(137.152);
}

// node_modules/d3-geo-projection/src/foucaut.js
function foucautRaw(lambda, phi) {
  var k2 = phi / 2, cosk = cos2(k2);
  return [2 * lambda / sqrtPi * cos2(phi) * cosk * cosk, sqrtPi * tan2(k2)];
}
foucautRaw.invert = function(x, y) {
  var k2 = atan3(y / sqrtPi), cosk = cos2(k2), phi = 2 * k2;
  return [x * sqrtPi / 2 / (cos2(phi) * cosk * cosk), phi];
};
function foucaut_default() {
  return projection(foucautRaw).scale(135.264);
}

// node_modules/d3-geo-projection/src/gilbert.js
function gilbertForward(point) {
  return [point[0] / 2, asin2(tan2(point[1] / 2 * radians2)) * degrees2];
}
function gilbertInvert(point) {
  return [point[0] * 2, 2 * atan3(sin2(point[1] * radians2)) * degrees2];
}
function gilbert_default(projectionType) {
  if (projectionType == null)
    projectionType = orthographic_default;
  var projection2 = projectionType(), equirectangular = equirectangular_default().scale(degrees2).precision(0).clipAngle(null).translate([0, 0]);
  function gilbert(point) {
    return projection2(gilbertForward(point));
  }
  if (projection2.invert)
    gilbert.invert = function(point) {
      return gilbertInvert(projection2.invert(point));
    };
  gilbert.stream = function(stream) {
    var s1 = projection2.stream(stream), s0 = equirectangular.stream({
      point: function(lambda, phi) {
        s1.point(lambda / 2, asin2(tan2(-phi / 2 * radians2)) * degrees2);
      },
      lineStart: function() {
        s1.lineStart();
      },
      lineEnd: function() {
        s1.lineEnd();
      },
      polygonStart: function() {
        s1.polygonStart();
      },
      polygonEnd: function() {
        s1.polygonEnd();
      }
    });
    s0.sphere = s1.sphere;
    return s0;
  };
  function property(name) {
    gilbert[name] = function() {
      return arguments.length ? (projection2[name].apply(projection2, arguments), gilbert) : projection2[name]();
    };
  }
  gilbert.rotate = function(_) {
    return arguments.length ? (equirectangular.rotate(_), gilbert) : equirectangular.rotate();
  };
  gilbert.center = function(_) {
    return arguments.length ? (projection2.center(gilbertForward(_)), gilbert) : gilbertInvert(projection2.center());
  };
  property("angle");
  property("clipAngle");
  property("clipExtent");
  property("fitExtent");
  property("fitHeight");
  property("fitSize");
  property("fitWidth");
  property("scale");
  property("translate");
  property("precision");
  return gilbert.scale(249.5);
}

// node_modules/d3-geo-projection/src/gingery.js
function gingeryRaw(rho, n) {
  var k2 = 2 * pi2 / n, rho2 = rho * rho;
  function forward(lambda, phi) {
    var p = azimuthalEquidistantRaw(lambda, phi), x = p[0], y = p[1], r2 = x * x + y * y;
    if (r2 > rho2) {
      var r = sqrt3(r2), theta = atan22(y, x), theta0 = k2 * round(theta / k2), alpha = theta - theta0, rhoCosAlpha = rho * cos2(alpha), k_ = (rho * sin2(alpha) - alpha * sin2(rhoCosAlpha)) / (halfPi2 - rhoCosAlpha), s_ = gingeryLength(alpha, k_), e = (pi2 - rho) / gingeryIntegrate(s_, rhoCosAlpha, pi2);
      x = r;
      var i = 50, delta;
      do {
        x -= delta = (rho + gingeryIntegrate(s_, rhoCosAlpha, x) * e - r) / (s_(x) * e);
      } while (abs2(delta) > epsilon3 && --i > 0);
      y = alpha * sin2(x);
      if (x < halfPi2)
        y -= k_ * (x - halfPi2);
      var s = sin2(theta0), c = cos2(theta0);
      p[0] = x * c - y * s;
      p[1] = x * s + y * c;
    }
    return p;
  }
  forward.invert = function(x, y) {
    var r2 = x * x + y * y;
    if (r2 > rho2) {
      var r = sqrt3(r2), theta = atan22(y, x), theta0 = k2 * round(theta / k2), dTheta = theta - theta0;
      x = r * cos2(dTheta);
      y = r * sin2(dTheta);
      var x_halfPi = x - halfPi2, sinx = sin2(x), alpha = y / sinx, delta = x < halfPi2 ? Infinity : 0, i = 10;
      while (true) {
        var rhosinAlpha = rho * sin2(alpha), rhoCosAlpha = rho * cos2(alpha), sinRhoCosAlpha = sin2(rhoCosAlpha), halfPi_RhoCosAlpha = halfPi2 - rhoCosAlpha, k_ = (rhosinAlpha - alpha * sinRhoCosAlpha) / halfPi_RhoCosAlpha, s_ = gingeryLength(alpha, k_);
        if (abs2(delta) < epsilon22 || !--i)
          break;
        alpha -= delta = (alpha * sinx - k_ * x_halfPi - y) / (sinx - x_halfPi * 2 * (halfPi_RhoCosAlpha * (rhoCosAlpha + alpha * rhosinAlpha * cos2(rhoCosAlpha) - sinRhoCosAlpha) - rhosinAlpha * (rhosinAlpha - alpha * sinRhoCosAlpha)) / (halfPi_RhoCosAlpha * halfPi_RhoCosAlpha));
      }
      r = rho + gingeryIntegrate(s_, rhoCosAlpha, x) * (pi2 - rho) / gingeryIntegrate(s_, rhoCosAlpha, pi2);
      theta = theta0 + alpha;
      x = r * cos2(theta);
      y = r * sin2(theta);
    }
    return azimuthalEquidistantRaw.invert(x, y);
  };
  return forward;
}
function gingeryLength(alpha, k2) {
  return function(x) {
    var y_ = alpha * cos2(x);
    if (x < halfPi2)
      y_ -= k2;
    return sqrt3(1 + y_ * y_);
  };
}
function gingeryIntegrate(f, a, b) {
  var n = 50, h = (b - a) / n, s = f(a) + f(b);
  for (var i = 1, x = a; i < n; ++i)
    s += 2 * f(x += h);
  return s * 0.5 * h;
}
function gingery_default() {
  var n = 6, rho = 30 * radians2, cRho = cos2(rho), sRho = sin2(rho), m = projectionMutator(gingeryRaw), p = m(rho, n), stream_ = p.stream, epsilon4 = 0.01, cr = -cos2(epsilon4 * radians2), sr = sin2(epsilon4 * radians2);
  p.radius = function(_) {
    if (!arguments.length)
      return rho * degrees2;
    cRho = cos2(rho = _ * radians2);
    sRho = sin2(rho);
    return m(rho, n);
  };
  p.lobes = function(_) {
    if (!arguments.length)
      return n;
    return m(rho, n = +_);
  };
  p.stream = function(stream) {
    var rotate = p.rotate(), rotateStream = stream_(stream), sphereStream = (p.rotate([0, 0]), stream_(stream));
    p.rotate(rotate);
    rotateStream.sphere = function() {
      sphereStream.polygonStart(), sphereStream.lineStart();
      for (var i = 0, delta = 2 * pi2 / n, phi = 0; i < n; ++i, phi -= delta) {
        sphereStream.point(atan22(sr * cos2(phi), cr) * degrees2, asin2(sr * sin2(phi)) * degrees2);
        sphereStream.point(atan22(sRho * cos2(phi - delta / 2), cRho) * degrees2, asin2(sRho * sin2(phi - delta / 2)) * degrees2);
      }
      sphereStream.lineEnd(), sphereStream.polygonEnd();
    };
    return rotateStream;
  };
  return p.rotate([90, -40]).scale(91.7095).clipAngle(180 - 1e-3);
}

// node_modules/d3-geo-projection/src/ginzburgPolyconic.js
function ginzburgPolyconic_default(a, b, c, d, e, f, g, h) {
  if (arguments.length < 8)
    h = 0;
  function forward(lambda, phi) {
    if (!phi)
      return [a * lambda / pi2, 0];
    var phi2 = phi * phi, xB = a + phi2 * (b + phi2 * (c + phi2 * d)), yB = phi * (e - 1 + phi2 * (f - h + phi2 * g)), m = (xB * xB + yB * yB) / (2 * yB), alpha = lambda * asin2(xB / m) / pi2;
    return [m * sin2(alpha), phi * (1 + phi2 * h) + m * (1 - cos2(alpha))];
  }
  forward.invert = function(x, y) {
    var lambda = pi2 * x / a, phi = y, deltaLambda, deltaPhi, i = 50;
    do {
      var phi2 = phi * phi, xB = a + phi2 * (b + phi2 * (c + phi2 * d)), yB = phi * (e - 1 + phi2 * (f - h + phi2 * g)), p = xB * xB + yB * yB, q = 2 * yB, m = p / q, m2 = m * m, dAlphadLambda = asin2(xB / m) / pi2, alpha = lambda * dAlphadLambda, xB2 = xB * xB, dxBdPhi = (2 * b + phi2 * (4 * c + phi2 * 6 * d)) * phi, dyBdPhi = e + phi2 * (3 * f + phi2 * 5 * g), dpdPhi = 2 * (xB * dxBdPhi + yB * (dyBdPhi - 1)), dqdPhi = 2 * (dyBdPhi - 1), dmdPhi = (dpdPhi * q - p * dqdPhi) / (q * q), cosAlpha = cos2(alpha), sinAlpha = sin2(alpha), mcosAlpha = m * cosAlpha, msinAlpha = m * sinAlpha, dAlphadPhi = lambda / pi2 * (1 / sqrt3(1 - xB2 / m2)) * (dxBdPhi * m - xB * dmdPhi) / m2, fx = msinAlpha - x, fy = phi * (1 + phi2 * h) + m - mcosAlpha - y, deltaxDeltaPhi = dmdPhi * sinAlpha + mcosAlpha * dAlphadPhi, deltaxDeltaLambda = mcosAlpha * dAlphadLambda, deltayDeltaPhi = 1 + dmdPhi - (dmdPhi * cosAlpha - msinAlpha * dAlphadPhi), deltayDeltaLambda = msinAlpha * dAlphadLambda, denominator = deltaxDeltaPhi * deltayDeltaLambda - deltayDeltaPhi * deltaxDeltaLambda;
      if (!denominator)
        break;
      lambda -= deltaLambda = (fy * deltaxDeltaPhi - fx * deltayDeltaPhi) / denominator;
      phi -= deltaPhi = (fx * deltayDeltaLambda - fy * deltaxDeltaLambda) / denominator;
    } while ((abs2(deltaLambda) > epsilon3 || abs2(deltaPhi) > epsilon3) && --i > 0);
    return [lambda, phi];
  };
  return forward;
}

// node_modules/d3-geo-projection/src/ginzburg4.js
var ginzburg4Raw = ginzburgPolyconic_default(2.8284, -1.6988, 0.75432, -0.18071, 1.76003, -0.38914, 0.042555);
function ginzburg4_default() {
  return projection(ginzburg4Raw).scale(149.995);
}

// node_modules/d3-geo-projection/src/ginzburg5.js
var ginzburg5Raw = ginzburgPolyconic_default(2.583819, -0.835827, 0.170354, -0.038094, 1.543313, -0.411435, 0.082742);
function ginzburg5_default() {
  return projection(ginzburg5Raw).scale(153.93);
}

// node_modules/d3-geo-projection/src/ginzburg6.js
var ginzburg6Raw = ginzburgPolyconic_default(5 / 6 * pi2, -0.62636, -0.0344, 0, 1.3493, -0.05524, 0, 0.045);
function ginzburg6_default() {
  return projection(ginzburg6Raw).scale(130.945);
}

// node_modules/d3-geo-projection/src/ginzburg8.js
function ginzburg8Raw(lambda, phi) {
  var lambda22 = lambda * lambda, phi2 = phi * phi;
  return [
    lambda * (1 - 0.162388 * phi2) * (0.87 - 952426e-9 * lambda22 * lambda22),
    phi * (1 + phi2 / 12)
  ];
}
ginzburg8Raw.invert = function(x, y) {
  var lambda = x, phi = y, i = 50, delta;
  do {
    var phi2 = phi * phi;
    phi -= delta = (phi * (1 + phi2 / 12) - y) / (1 + phi2 / 4);
  } while (abs2(delta) > epsilon3 && --i > 0);
  i = 50;
  x /= 1 - 0.162388 * phi2;
  do {
    var lambda4 = (lambda4 = lambda * lambda) * lambda4;
    lambda -= delta = (lambda * (0.87 - 952426e-9 * lambda4) - x) / (0.87 - 476213e-8 * lambda4);
  } while (abs2(delta) > epsilon3 && --i > 0);
  return [lambda, phi];
};
function ginzburg8_default() {
  return projection(ginzburg8Raw).scale(131.747);
}

// node_modules/d3-geo-projection/src/ginzburg9.js
var ginzburg9Raw = ginzburgPolyconic_default(2.6516, -0.76534, 0.19123, -0.047094, 1.36289, -0.13965, 0.031762);
function ginzburg9_default() {
  return projection(ginzburg9Raw).scale(131.087);
}

// node_modules/d3-geo-projection/src/square.js
function square_default(project) {
  var dx = project(halfPi2, 0)[0] - project(-halfPi2, 0)[0];
  function projectSquare(lambda, phi) {
    var s = lambda > 0 ? -0.5 : 0.5, point = project(lambda + s * pi2, phi);
    point[0] -= s * dx;
    return point;
  }
  if (project.invert)
    projectSquare.invert = function(x, y) {
      var s = x > 0 ? -0.5 : 0.5, location = project.invert(x + s * dx, y), lambda = location[0] - s * pi2;
      if (lambda < -pi2)
        lambda += 2 * pi2;
      else if (lambda > pi2)
        lambda -= 2 * pi2;
      location[0] = lambda;
      return location;
    };
  return projectSquare;
}

// node_modules/d3-geo-projection/src/gringorten.js
function gringortenRaw(lambda, phi) {
  var sLambda = sign2(lambda), sPhi = sign2(phi), cosPhi = cos2(phi), x = cos2(lambda) * cosPhi, y = sin2(lambda) * cosPhi, z = sin2(sPhi * phi);
  lambda = abs2(atan22(y, z));
  phi = asin2(x);
  if (abs2(lambda - halfPi2) > epsilon3)
    lambda %= halfPi2;
  var point = gringortenHexadecant(lambda > pi2 / 4 ? halfPi2 - lambda : lambda, phi);
  if (lambda > pi2 / 4)
    z = point[0], point[0] = -point[1], point[1] = -z;
  return point[0] *= sLambda, point[1] *= -sPhi, point;
}
gringortenRaw.invert = function(x, y) {
  if (abs2(x) > 1)
    x = sign2(x) * 2 - x;
  if (abs2(y) > 1)
    y = sign2(y) * 2 - y;
  var sx = sign2(x), sy = sign2(y), x03 = -sx * x, y03 = -sy * y, t = y03 / x03 < 1, p = gringortenHexadecantInvert(t ? y03 : x03, t ? x03 : y03), lambda = p[0], phi = p[1], cosPhi = cos2(phi);
  if (t)
    lambda = -halfPi2 - lambda;
  return [sx * (atan22(sin2(lambda) * cosPhi, -sin2(phi)) + pi2), sy * asin2(cos2(lambda) * cosPhi)];
};
function gringortenHexadecant(lambda, phi) {
  if (phi === halfPi2)
    return [0, 0];
  var sinPhi = sin2(phi), r = sinPhi * sinPhi, r2 = r * r, j = 1 + r2, k2 = 1 + 3 * r2, q = 1 - r2, z = asin2(1 / sqrt3(j)), v = q + r * j * z, p2 = (1 - sinPhi) / v, p = sqrt3(p2), a2 = p2 * j, a = sqrt3(a2), h = p * q, x, i;
  if (lambda === 0)
    return [0, -(h + r * a)];
  var cosPhi = cos2(phi), secPhi = 1 / cosPhi, drdPhi = 2 * sinPhi * cosPhi, dvdPhi = (-3 * r + z * k2) * drdPhi, dp2dPhi = (-v * cosPhi - (1 - sinPhi) * dvdPhi) / (v * v), dpdPhi = 0.5 * dp2dPhi / p, dhdPhi = q * dpdPhi - 2 * r * p * drdPhi, dra2dPhi = r * j * dp2dPhi + p2 * k2 * drdPhi, mu = -secPhi * drdPhi, nu = -secPhi * dra2dPhi, zeta = -2 * secPhi * dhdPhi, lambda12 = 4 * lambda / pi2, delta;
  if (lambda > 0.222 * pi2 || phi < pi2 / 4 && lambda > 0.175 * pi2) {
    x = (h + r * sqrt3(a2 * (1 + r2) - h * h)) / (1 + r2);
    if (lambda > pi2 / 4)
      return [x, x];
    var x12 = x, x03 = 0.5 * x;
    x = 0.5 * (x03 + x12), i = 50;
    do {
      var g = sqrt3(a2 - x * x), f = x * (zeta + mu * g) + nu * asin2(x / a) - lambda12;
      if (!f)
        break;
      if (f < 0)
        x03 = x;
      else
        x12 = x;
      x = 0.5 * (x03 + x12);
    } while (abs2(x12 - x03) > epsilon3 && --i > 0);
  } else {
    x = epsilon3, i = 25;
    do {
      var x2 = x * x, g2 = sqrt3(a2 - x2), zetaMug = zeta + mu * g2, f2 = x * zetaMug + nu * asin2(x / a) - lambda12, df = zetaMug + (nu - mu * x2) / g2;
      x -= delta = g2 ? f2 / df : 0;
    } while (abs2(delta) > epsilon3 && --i > 0);
  }
  return [x, -h - r * sqrt3(a2 - x * x)];
}
function gringortenHexadecantInvert(x, y) {
  var x03 = 0, x12 = 1, r = 0.5, i = 50;
  while (true) {
    var r2 = r * r, sinPhi = sqrt3(r), z = asin2(1 / sqrt3(1 + r2)), v = 1 - r2 + r * (1 + r2) * z, p2 = (1 - sinPhi) / v, p = sqrt3(p2), a2 = p2 * (1 + r2), h = p * (1 - r2), g2 = a2 - x * x, g = sqrt3(g2), y03 = y + h + r * g;
    if (abs2(x12 - x03) < epsilon22 || --i === 0 || y03 === 0)
      break;
    if (y03 > 0)
      x03 = r;
    else
      x12 = r;
    r = 0.5 * (x03 + x12);
  }
  if (!i)
    return null;
  var phi = asin2(sinPhi), cosPhi = cos2(phi), secPhi = 1 / cosPhi, drdPhi = 2 * sinPhi * cosPhi, dvdPhi = (-3 * r + z * (1 + 3 * r2)) * drdPhi, dp2dPhi = (-v * cosPhi - (1 - sinPhi) * dvdPhi) / (v * v), dpdPhi = 0.5 * dp2dPhi / p, dhdPhi = (1 - r2) * dpdPhi - 2 * r * p * drdPhi, zeta = -2 * secPhi * dhdPhi, mu = -secPhi * drdPhi, nu = -secPhi * (r * (1 + r2) * dp2dPhi + p2 * (1 + 3 * r2) * drdPhi);
  return [pi2 / 4 * (x * (zeta + mu * g) + nu * asin2(x / sqrt3(a2))), phi];
}
function gringorten_default() {
  return projection(square_default(gringortenRaw)).scale(239.75);
}

// node_modules/d3-geo-projection/src/elliptic.js
function ellipticJi(u, v, m) {
  var a, b, c;
  if (!u) {
    b = ellipticJ(v, 1 - m);
    return [
      [0, b[0] / b[1]],
      [1 / b[1], 0],
      [b[2] / b[1], 0]
    ];
  }
  a = ellipticJ(u, m);
  if (!v)
    return [[a[0], 0], [a[1], 0], [a[2], 0]];
  b = ellipticJ(v, 1 - m);
  c = b[1] * b[1] + m * a[0] * a[0] * b[0] * b[0];
  return [
    [a[0] * b[2] / c, a[1] * a[2] * b[0] * b[1] / c],
    [a[1] * b[1] / c, -a[0] * a[2] * b[0] * b[2] / c],
    [a[2] * b[1] * b[2] / c, -m * a[0] * a[1] * b[0] / c]
  ];
}
function ellipticJ(u, m) {
  var ai, b, phi, t, twon;
  if (m < epsilon3) {
    t = sin2(u);
    b = cos2(u);
    ai = m * (u - t * b) / 4;
    return [
      t - ai * b,
      b + ai * t,
      1 - m * t * t / 2,
      u - ai
    ];
  }
  if (m >= 1 - epsilon3) {
    ai = (1 - m) / 4;
    b = cosh(u);
    t = tanh(u);
    phi = 1 / b;
    twon = b * sinh(u);
    return [
      t + ai * (twon - u) / (b * b),
      phi - ai * t * phi * (twon - u),
      phi + ai * t * phi * (twon + u),
      2 * atan3(exp2(u)) - halfPi2 + ai * (twon - u) / b
    ];
  }
  var a = [1, 0, 0, 0, 0, 0, 0, 0, 0], c = [sqrt3(m), 0, 0, 0, 0, 0, 0, 0, 0], i = 0;
  b = sqrt3(1 - m);
  twon = 1;
  while (abs2(c[i] / a[i]) > epsilon3 && i < 8) {
    ai = a[i++];
    c[i] = (ai - b) / 2;
    a[i] = (ai + b) / 2;
    b = sqrt3(ai * b);
    twon *= 2;
  }
  phi = twon * a[i] * u;
  do {
    t = c[i] * sin2(b = phi) / a[i];
    phi = (asin2(t) + phi) / 2;
  } while (--i);
  return [sin2(phi), t = cos2(phi), t / cos2(phi - b), phi];
}
function ellipticFi(phi, psi, m) {
  var r = abs2(phi), i = abs2(psi), sinhPsi = sinh(i);
  if (r) {
    var cscPhi = 1 / sin2(r), cotPhi2 = 1 / (tan2(r) * tan2(r)), b = -(cotPhi2 + m * (sinhPsi * sinhPsi * cscPhi * cscPhi) - 1 + m), c = (m - 1) * cotPhi2, cotLambda2 = (-b + sqrt3(b * b - 4 * c)) / 2;
    return [
      ellipticF(atan3(1 / sqrt3(cotLambda2)), m) * sign2(phi),
      ellipticF(atan3(sqrt3((cotLambda2 / cotPhi2 - 1) / m)), 1 - m) * sign2(psi)
    ];
  }
  return [
    0,
    ellipticF(atan3(sinhPsi), 1 - m) * sign2(psi)
  ];
}
function ellipticF(phi, m) {
  if (!m)
    return phi;
  if (m === 1)
    return log2(tan2(phi / 2 + quarterPi2));
  var a = 1, b = sqrt3(1 - m), c = sqrt3(m);
  for (var i = 0; abs2(c) > epsilon3; i++) {
    if (phi % pi2) {
      var dPhi = atan3(b * tan2(phi) / a);
      if (dPhi < 0)
        dPhi += pi2;
      phi += dPhi + ~~(phi / pi2) * pi2;
    } else
      phi += phi;
    c = (a + b) / 2;
    b = sqrt3(a * b);
    c = ((a = c) - b) / 2;
  }
  return phi / (pow2(2, i) * a);
}

// node_modules/d3-geo-projection/src/guyou.js
function guyouRaw(lambda, phi) {
  var k_ = (sqrt2 - 1) / (sqrt2 + 1), k2 = sqrt3(1 - k_ * k_), K3 = ellipticF(halfPi2, k2 * k2), f = -1, psi = log2(tan2(pi2 / 4 + abs2(phi) / 2)), r = exp2(f * psi) / sqrt3(k_), at = guyouComplexAtan(r * cos2(f * lambda), r * sin2(f * lambda)), t = ellipticFi(at[0], at[1], k2 * k2);
  return [-t[1], (phi >= 0 ? 1 : -1) * (0.5 * K3 - t[0])];
}
function guyouComplexAtan(x, y) {
  var x2 = x * x, y_1 = y + 1, t = 1 - x2 - y * y;
  return [
    0.5 * ((x >= 0 ? halfPi2 : -halfPi2) - atan22(t, 2 * x)),
    -0.25 * log2(t * t + 4 * x2) + 0.5 * log2(y_1 * y_1 + x2)
  ];
}
function guyouComplexDivide(a, b) {
  var denominator = b[0] * b[0] + b[1] * b[1];
  return [
    (a[0] * b[0] + a[1] * b[1]) / denominator,
    (a[1] * b[0] - a[0] * b[1]) / denominator
  ];
}
guyouRaw.invert = function(x, y) {
  var k_ = (sqrt2 - 1) / (sqrt2 + 1), k2 = sqrt3(1 - k_ * k_), K3 = ellipticF(halfPi2, k2 * k2), f = -1, j = ellipticJi(0.5 * K3 - y, -x, k2 * k2), tn = guyouComplexDivide(j[0], j[1]), lambda = atan22(tn[1], tn[0]) / f;
  return [
    lambda,
    2 * atan3(exp2(0.5 / f * log2(k_ * tn[0] * tn[0] + k_ * tn[1] * tn[1]))) - halfPi2
  ];
};
function guyou_default() {
  return projection(square_default(guyouRaw)).scale(151.496);
}

// node_modules/d3-geo-projection/src/hammerRetroazimuthal.js
function hammerRetroazimuthalRaw(phi03) {
  var sinPhi02 = sin2(phi03), cosPhi02 = cos2(phi03), rotate = hammerRetroazimuthalRotation(phi03);
  rotate.invert = hammerRetroazimuthalRotation(-phi03);
  function forward(lambda, phi) {
    var p = rotate(lambda, phi);
    lambda = p[0], phi = p[1];
    var sinPhi = sin2(phi), cosPhi = cos2(phi), cosLambda = cos2(lambda), z = acos2(sinPhi02 * sinPhi + cosPhi02 * cosPhi * cosLambda), sinz = sin2(z), K3 = abs2(sinz) > epsilon3 ? z / sinz : 1;
    return [
      K3 * cosPhi02 * sin2(lambda),
      (abs2(lambda) > halfPi2 ? K3 : -K3) * (sinPhi02 * cosPhi - cosPhi02 * sinPhi * cosLambda)
    ];
  }
  forward.invert = function(x, y) {
    var rho = sqrt3(x * x + y * y), sinz = -sin2(rho), cosz = cos2(rho), a = rho * cosz, b = -y * sinz, c = rho * sinPhi02, d = sqrt3(a * a + b * b - c * c), phi = atan22(a * c + b * d, b * c - a * d), lambda = (rho > halfPi2 ? -1 : 1) * atan22(x * sinz, rho * cos2(phi) * cosz + y * sin2(phi) * sinz);
    return rotate.invert(lambda, phi);
  };
  return forward;
}
function hammerRetroazimuthalRotation(phi03) {
  var sinPhi02 = sin2(phi03), cosPhi02 = cos2(phi03);
  return function(lambda, phi) {
    var cosPhi = cos2(phi), x = cos2(lambda) * cosPhi, y = sin2(lambda) * cosPhi, z = sin2(phi);
    return [
      atan22(y, x * cosPhi02 - z * sinPhi02),
      asin2(z * cosPhi02 + x * sinPhi02)
    ];
  };
}
function hammerRetroazimuthal_default() {
  var phi03 = 0, m = projectionMutator(hammerRetroazimuthalRaw), p = m(phi03), rotate_ = p.rotate, stream_ = p.stream, circle = circle_default();
  p.parallel = function(_) {
    if (!arguments.length)
      return phi03 * degrees2;
    var r = p.rotate();
    return m(phi03 = _ * radians2).rotate(r);
  };
  p.rotate = function(_) {
    if (!arguments.length)
      return _ = rotate_.call(p), _[1] += phi03 * degrees2, _;
    rotate_.call(p, [_[0], _[1] - phi03 * degrees2]);
    circle.center([-_[0], -_[1]]);
    return p;
  };
  p.stream = function(stream) {
    stream = stream_(stream);
    stream.sphere = function() {
      stream.polygonStart();
      var epsilon4 = 0.01, ring = circle.radius(90 - epsilon4)().coordinates[0], n = ring.length - 1, i = -1, p2;
      stream.lineStart();
      while (++i < n)
        stream.point((p2 = ring[i])[0], p2[1]);
      stream.lineEnd();
      ring = circle.radius(90 + epsilon4)().coordinates[0];
      n = ring.length - 1;
      stream.lineStart();
      while (--i >= 0)
        stream.point((p2 = ring[i])[0], p2[1]);
      stream.lineEnd();
      stream.polygonEnd();
    };
    return stream;
  };
  return p.scale(79.4187).parallel(45).clipAngle(180 - 1e-3);
}

// node_modules/d3-geo-projection/src/healpix.js
var K = 3;
var healpixParallel = asin2(1 - 1 / K) * degrees2;
var healpixLambert = cylindricalEqualAreaRaw2(0);
function healpixRaw(H) {
  var phi03 = healpixParallel * radians2, dx = collignonRaw(pi2, phi03)[0] - collignonRaw(-pi2, phi03)[0], y03 = healpixLambert(0, phi03)[1], y12 = collignonRaw(0, phi03)[1], dy1 = sqrtPi - y12, k2 = tau2 / H, w2 = 4 / tau2, h = y03 + dy1 * dy1 * 4 / tau2;
  function forward(lambda, phi) {
    var point, phi2 = abs2(phi);
    if (phi2 > phi03) {
      var i = min(H - 1, max(0, floor((lambda + pi2) / k2)));
      lambda += pi2 * (H - 1) / H - i * k2;
      point = collignonRaw(lambda, phi2);
      point[0] = point[0] * tau2 / dx - tau2 * (H - 1) / (2 * H) + i * tau2 / H;
      point[1] = y03 + (point[1] - y12) * 4 * dy1 / tau2;
      if (phi < 0)
        point[1] = -point[1];
    } else {
      point = healpixLambert(lambda, phi);
    }
    point[0] *= w2, point[1] /= h;
    return point;
  }
  forward.invert = function(x, y) {
    x /= w2, y *= h;
    var y2 = abs2(y);
    if (y2 > y03) {
      var i = min(H - 1, max(0, floor((x + pi2) / k2)));
      x = (x + pi2 * (H - 1) / H - i * k2) * dx / tau2;
      var point = collignonRaw.invert(x, 0.25 * (y2 - y03) * tau2 / dy1 + y12);
      point[0] -= pi2 * (H - 1) / H - i * k2;
      if (y < 0)
        point[1] = -point[1];
      return point;
    }
    return healpixLambert.invert(x, y);
  };
  return forward;
}
function sphereTop(x, i) {
  return [x, i & 1 ? 90 - epsilon3 : healpixParallel];
}
function sphereBottom(x, i) {
  return [x, i & 1 ? -90 + epsilon3 : -healpixParallel];
}
function sphereNudge(d) {
  return [d[0] * (1 - epsilon3), d[1]];
}
function sphere(step) {
  var c = [].concat(
    range(-180, 180 + step / 2, step).map(sphereTop),
    range(180, -180 - step / 2, -step).map(sphereBottom)
  );
  return {
    type: "Polygon",
    coordinates: [step === 180 ? c.map(sphereNudge) : c]
  };
}
function healpix_default() {
  var H = 4, m = projectionMutator(healpixRaw), p = m(H), stream_ = p.stream;
  p.lobes = function(_) {
    return arguments.length ? m(H = +_) : H;
  };
  p.stream = function(stream) {
    var rotate = p.rotate(), rotateStream = stream_(stream), sphereStream = (p.rotate([0, 0]), stream_(stream));
    p.rotate(rotate);
    rotateStream.sphere = function() {
      stream_default(sphere(180 / H), sphereStream);
    };
    return rotateStream;
  };
  return p.scale(239.75);
}

// node_modules/d3-geo-projection/src/hill.js
function hillRaw(K3) {
  var L = 1 + K3, sinBt = sin2(1 / L), Bt = asin2(sinBt), A5 = 2 * sqrt3(pi2 / (B2 = pi2 + 4 * Bt * L)), B2, rho0 = 0.5 * A5 * (L + sqrt3(K3 * (2 + K3))), K22 = K3 * K3, L2 = L * L;
  function forward(lambda, phi) {
    var t = 1 - sin2(phi), rho, omega;
    if (t && t < 2) {
      var theta = halfPi2 - phi, i = 25, delta;
      do {
        var sinTheta = sin2(theta), cosTheta = cos2(theta), Bt_Bt1 = Bt + atan22(sinTheta, L - cosTheta), C = 1 + L2 - 2 * L * cosTheta;
        theta -= delta = (theta - K22 * Bt - L * sinTheta + C * Bt_Bt1 - 0.5 * t * B2) / (2 * L * sinTheta * Bt_Bt1);
      } while (abs2(delta) > epsilon22 && --i > 0);
      rho = A5 * sqrt3(C);
      omega = lambda * Bt_Bt1 / pi2;
    } else {
      rho = A5 * (K3 + t);
      omega = lambda * Bt / pi2;
    }
    return [
      rho * sin2(omega),
      rho0 - rho * cos2(omega)
    ];
  }
  forward.invert = function(x, y) {
    var rho2 = x * x + (y -= rho0) * y, cosTheta = (1 + L2 - rho2 / (A5 * A5)) / (2 * L), theta = acos2(cosTheta), sinTheta = sin2(theta), Bt_Bt1 = Bt + atan22(sinTheta, L - cosTheta);
    return [
      asin2(x / sqrt3(rho2)) * pi2 / Bt_Bt1,
      asin2(1 - 2 * (theta - K22 * Bt - L * sinTheta + (1 + L2 - 2 * L * cosTheta) * Bt_Bt1) / B2)
    ];
  };
  return forward;
}
function hill_default() {
  var K3 = 1, m = projectionMutator(hillRaw), p = m(K3);
  p.ratio = function(_) {
    return arguments.length ? m(K3 = +_) : K3;
  };
  return p.scale(167.774).center([0, 18.67]);
}

// node_modules/d3-geo-projection/src/sinuMollweide.js
var sinuMollweidePhi = 0.7109889596207567;
var sinuMollweideY = 0.0528035274542;
function sinuMollweideRaw(lambda, phi) {
  return phi > -sinuMollweidePhi ? (lambda = mollweideRaw(lambda, phi), lambda[1] += sinuMollweideY, lambda) : sinusoidalRaw(lambda, phi);
}
sinuMollweideRaw.invert = function(x, y) {
  return y > -sinuMollweidePhi ? mollweideRaw.invert(x, y - sinuMollweideY) : sinusoidalRaw.invert(x, y);
};
function sinuMollweide_default() {
  return projection(sinuMollweideRaw).rotate([-20, -55]).scale(164.263).center([0, -5.4036]);
}

// node_modules/d3-geo-projection/src/homolosine.js
function homolosineRaw(lambda, phi) {
  return abs2(phi) > sinuMollweidePhi ? (lambda = mollweideRaw(lambda, phi), lambda[1] -= phi > 0 ? sinuMollweideY : -sinuMollweideY, lambda) : sinusoidalRaw(lambda, phi);
}
homolosineRaw.invert = function(x, y) {
  return abs2(y) > sinuMollweidePhi ? mollweideRaw.invert(x, y + (y > 0 ? sinuMollweideY : -sinuMollweideY)) : sinusoidalRaw.invert(x, y);
};
function homolosine_default() {
  return projection(homolosineRaw).scale(152.63);
}

// node_modules/d3-geo-projection/src/interrupted/index.js
function pointEqual(a, b) {
  return abs2(a[0] - b[0]) < epsilon3 && abs2(a[1] - b[1]) < epsilon3;
}
function interpolateLine(coordinates, m) {
  var i = -1, n = coordinates.length, p02 = coordinates[0], p1, dx, dy, resampled = [];
  while (++i < n) {
    p1 = coordinates[i];
    dx = (p1[0] - p02[0]) / m;
    dy = (p1[1] - p02[1]) / m;
    for (var j = 0; j < m; ++j)
      resampled.push([p02[0] + j * dx, p02[1] + j * dy]);
    p02 = p1;
  }
  resampled.push(p1);
  return resampled;
}
function interpolateSphere(lobes7) {
  var coordinates = [], lobe, lambda03, phi03, phi12, lambda22, phi2, i, n = lobes7[0].length;
  for (i = 0; i < n; ++i) {
    lobe = lobes7[0][i];
    lambda03 = lobe[0][0], phi03 = lobe[0][1], phi12 = lobe[1][1];
    lambda22 = lobe[2][0], phi2 = lobe[2][1];
    coordinates.push(interpolateLine([
      [lambda03 + epsilon3, phi03 + epsilon3],
      [lambda03 + epsilon3, phi12 - epsilon3],
      [lambda22 - epsilon3, phi12 - epsilon3],
      [lambda22 - epsilon3, phi2 + epsilon3]
    ], 30));
  }
  for (i = lobes7[1].length - 1; i >= 0; --i) {
    lobe = lobes7[1][i];
    lambda03 = lobe[0][0], phi03 = lobe[0][1], phi12 = lobe[1][1];
    lambda22 = lobe[2][0], phi2 = lobe[2][1];
    coordinates.push(interpolateLine([
      [lambda22 - epsilon3, phi2 - epsilon3],
      [lambda22 - epsilon3, phi12 + epsilon3],
      [lambda03 + epsilon3, phi12 + epsilon3],
      [lambda03 + epsilon3, phi03 - epsilon3]
    ], 30));
  }
  return {
    type: "Polygon",
    coordinates: [merge(coordinates)]
  };
}
function interrupted_default(project, lobes7, inverse2) {
  var sphere2, bounds;
  function forward(lambda, phi) {
    var sign3 = phi < 0 ? -1 : 1, lobe = lobes7[+(phi < 0)];
    for (var i = 0, n = lobe.length - 1; i < n && lambda > lobe[i][2][0]; ++i)
      ;
    var p2 = project(lambda - lobe[i][1][0], phi);
    p2[0] += project(lobe[i][1][0], sign3 * phi > sign3 * lobe[i][0][1] ? lobe[i][0][1] : phi)[0];
    return p2;
  }
  if (inverse2) {
    forward.invert = inverse2(forward);
  } else if (project.invert) {
    forward.invert = function(x, y) {
      var bound = bounds[+(y < 0)], lobe = lobes7[+(y < 0)];
      for (var i = 0, n = bound.length; i < n; ++i) {
        var b = bound[i];
        if (b[0][0] <= x && x < b[1][0] && b[0][1] <= y && y < b[1][1]) {
          var p2 = project.invert(x - project(lobe[i][1][0], 0)[0], y);
          p2[0] += lobe[i][1][0];
          return pointEqual(forward(p2[0], p2[1]), [x, y]) ? p2 : null;
        }
      }
    };
  }
  var p = projection(forward), stream_ = p.stream;
  p.stream = function(stream) {
    var rotate = p.rotate(), rotateStream = stream_(stream), sphereStream = (p.rotate([0, 0]), stream_(stream));
    p.rotate(rotate);
    rotateStream.sphere = function() {
      stream_default(sphere2, sphereStream);
    };
    return rotateStream;
  };
  p.lobes = function(_) {
    if (!arguments.length)
      return lobes7.map(function(lobe) {
        return lobe.map(function(l) {
          return [
            [l[0][0] * degrees2, l[0][1] * degrees2],
            [l[1][0] * degrees2, l[1][1] * degrees2],
            [l[2][0] * degrees2, l[2][1] * degrees2]
          ];
        });
      });
    sphere2 = interpolateSphere(_);
    lobes7 = _.map(function(lobe) {
      return lobe.map(function(l) {
        return [
          [l[0][0] * radians2, l[0][1] * radians2],
          [l[1][0] * radians2, l[1][1] * radians2],
          [l[2][0] * radians2, l[2][1] * radians2]
        ];
      });
    });
    bounds = lobes7.map(function(lobe) {
      return lobe.map(function(l) {
        var x03 = project(l[0][0], l[0][1])[0], x12 = project(l[2][0], l[2][1])[0], y03 = project(l[1][0], l[0][1])[1], y12 = project(l[1][0], l[1][1])[1], t;
        if (y03 > y12)
          t = y03, y03 = y12, y12 = t;
        return [[x03, y03], [x12, y12]];
      });
    });
    return p;
  };
  if (lobes7 != null)
    p.lobes(lobes7);
  return p;
}

// node_modules/d3-geo-projection/src/interrupted/boggs.js
var lobes = [[
  // northern hemisphere
  [[-180, 0], [-100, 90], [-40, 0]],
  [[-40, 0], [30, 90], [180, 0]]
], [
  // southern hemisphere
  [[-180, 0], [-160, -90], [-100, 0]],
  [[-100, 0], [-60, -90], [-20, 0]],
  [[-20, 0], [20, -90], [80, 0]],
  [[80, 0], [140, -90], [180, 0]]
]];
function boggs_default2() {
  return interrupted_default(boggsRaw, lobes).scale(160.857);
}

// node_modules/d3-geo-projection/src/interrupted/homolosine.js
var lobes2 = [[
  // northern hemisphere
  [[-180, 0], [-100, 90], [-40, 0]],
  [[-40, 0], [30, 90], [180, 0]]
], [
  // southern hemisphere
  [[-180, 0], [-160, -90], [-100, 0]],
  [[-100, 0], [-60, -90], [-20, 0]],
  [[-20, 0], [20, -90], [80, 0]],
  [[80, 0], [140, -90], [180, 0]]
]];
function homolosine_default2() {
  return interrupted_default(homolosineRaw, lobes2).scale(152.63);
}

// node_modules/d3-geo-projection/src/interrupted/mollweide.js
var lobes3 = [[
  // northern hemisphere
  [[-180, 0], [-100, 90], [-40, 0]],
  [[-40, 0], [30, 90], [180, 0]]
], [
  // southern hemisphere
  [[-180, 0], [-160, -90], [-100, 0]],
  [[-100, 0], [-60, -90], [-20, 0]],
  [[-20, 0], [20, -90], [80, 0]],
  [[80, 0], [140, -90], [180, 0]]
]];
function mollweide_default2() {
  return interrupted_default(mollweideRaw, lobes3).scale(169.529);
}

// node_modules/d3-geo-projection/src/interrupted/mollweideHemispheres.js
var lobes4 = [[
  // northern hemisphere
  [[-180, 0], [-90, 90], [0, 0]],
  [[0, 0], [90, 90], [180, 0]]
], [
  // southern hemisphere
  [[-180, 0], [-90, -90], [0, 0]],
  [[0, 0], [90, -90], [180, 0]]
]];
function mollweideHemispheres_default() {
  return interrupted_default(mollweideRaw, lobes4).scale(169.529).rotate([20, 0]);
}

// node_modules/d3-geo-projection/src/interrupted/sinuMollweide.js
var lobes5 = [[
  // northern hemisphere
  [[-180, 35], [-30, 90], [0, 35]],
  [[0, 35], [30, 90], [180, 35]]
], [
  // southern hemisphere
  [[-180, -10], [-102, -90], [-65, -10]],
  [[-65, -10], [5, -90], [77, -10]],
  [[77, -10], [103, -90], [180, -10]]
]];
function sinuMollweide_default2() {
  return interrupted_default(sinuMollweideRaw, lobes5, solve2d).rotate([-20, -55]).scale(164.263).center([0, -5.4036]);
}

// node_modules/d3-geo-projection/src/interrupted/sinusoidal.js
var lobes6 = [[
  // northern hemisphere
  [[-180, 0], [-110, 90], [-40, 0]],
  [[-40, 0], [0, 90], [40, 0]],
  [[40, 0], [110, 90], [180, 0]]
], [
  // southern hemisphere
  [[-180, 0], [-110, -90], [-40, 0]],
  [[-40, 0], [0, -90], [40, 0]],
  [[40, 0], [110, -90], [180, 0]]
]];
function sinusoidal_default2() {
  return interrupted_default(sinusoidalRaw, lobes6).scale(152.63).rotate([-20, 0]);
}

// node_modules/d3-geo-projection/src/kavrayskiy7.js
function kavrayskiy7Raw(lambda, phi) {
  return [3 / tau2 * lambda * sqrt3(pi2 * pi2 / 3 - phi * phi), phi];
}
kavrayskiy7Raw.invert = function(x, y) {
  return [tau2 / 3 * x / sqrt3(pi2 * pi2 / 3 - y * y), y];
};
function kavrayskiy7_default() {
  return projection(kavrayskiy7Raw).scale(158.837);
}

// node_modules/d3-geo-projection/src/lagrange.js
function lagrangeRaw(n) {
  function forward(lambda, phi) {
    if (abs2(abs2(phi) - halfPi2) < epsilon3)
      return [0, phi < 0 ? -2 : 2];
    var sinPhi = sin2(phi), v = pow2((1 + sinPhi) / (1 - sinPhi), n / 2), c = 0.5 * (v + 1 / v) + cos2(lambda *= n);
    return [
      2 * sin2(lambda) / c,
      (v - 1 / v) / c
    ];
  }
  forward.invert = function(x, y) {
    var y03 = abs2(y);
    if (abs2(y03 - 2) < epsilon3)
      return x ? null : [0, sign2(y) * halfPi2];
    if (y03 > 2)
      return null;
    x /= 2, y /= 2;
    var x2 = x * x, y2 = y * y, t = 2 * y / (1 + x2 + y2);
    t = pow2((1 + t) / (1 - t), 1 / n);
    return [
      atan22(2 * x, 1 - x2 - y2) / n,
      asin2((t - 1) / (t + 1))
    ];
  };
  return forward;
}
function lagrange_default() {
  var n = 0.5, m = projectionMutator(lagrangeRaw), p = m(n);
  p.spacing = function(_) {
    return arguments.length ? m(n = +_) : n;
  };
  return p.scale(124.75);
}

// node_modules/d3-geo-projection/src/larrivee.js
var pi_sqrt2 = pi2 / sqrt2;
function larriveeRaw(lambda, phi) {
  return [
    lambda * (1 + sqrt3(cos2(phi))) / 2,
    phi / (cos2(phi / 2) * cos2(lambda / 6))
  ];
}
larriveeRaw.invert = function(x, y) {
  var x03 = abs2(x), y03 = abs2(y), lambda = epsilon3, phi = halfPi2;
  if (y03 < pi_sqrt2)
    phi *= y03 / pi_sqrt2;
  else
    lambda += 6 * acos2(pi_sqrt2 / y03);
  for (var i = 0; i < 25; i++) {
    var sinPhi = sin2(phi), sqrtcosPhi = sqrt3(cos2(phi)), sinPhi_2 = sin2(phi / 2), cosPhi_2 = cos2(phi / 2), sinLambda_6 = sin2(lambda / 6), cosLambda_6 = cos2(lambda / 6), f0 = 0.5 * lambda * (1 + sqrtcosPhi) - x03, f1 = phi / (cosPhi_2 * cosLambda_6) - y03, df0dPhi = sqrtcosPhi ? -0.25 * lambda * sinPhi / sqrtcosPhi : 0, df0dLambda = 0.5 * (1 + sqrtcosPhi), df1dPhi = (1 + 0.5 * phi * sinPhi_2 / cosPhi_2) / (cosPhi_2 * cosLambda_6), df1dLambda = phi / cosPhi_2 * (sinLambda_6 / 6) / (cosLambda_6 * cosLambda_6), denom = df0dPhi * df1dLambda - df1dPhi * df0dLambda, dPhi = (f0 * df1dLambda - f1 * df0dLambda) / denom, dLambda = (f1 * df0dPhi - f0 * df1dPhi) / denom;
    phi -= dPhi;
    lambda -= dLambda;
    if (abs2(dPhi) < epsilon3 && abs2(dLambda) < epsilon3)
      break;
  }
  return [x < 0 ? -lambda : lambda, y < 0 ? -phi : phi];
};
function larrivee_default() {
  return projection(larriveeRaw).scale(97.2672);
}

// node_modules/d3-geo-projection/src/laskowski.js
function laskowskiRaw(lambda, phi) {
  var lambda22 = lambda * lambda, phi2 = phi * phi;
  return [
    lambda * (0.975534 + phi2 * (-0.119161 + lambda22 * -0.0143059 + phi2 * -0.0547009)),
    phi * (1.00384 + lambda22 * (0.0802894 + phi2 * -0.02855 + lambda22 * 199025e-9) + phi2 * (0.0998909 + phi2 * -0.0491032))
  ];
}
laskowskiRaw.invert = function(x, y) {
  var lambda = sign2(x) * pi2, phi = y / 2, i = 50;
  do {
    var lambda22 = lambda * lambda, phi2 = phi * phi, lambdaPhi = lambda * phi, fx = lambda * (0.975534 + phi2 * (-0.119161 + lambda22 * -0.0143059 + phi2 * -0.0547009)) - x, fy = phi * (1.00384 + lambda22 * (0.0802894 + phi2 * -0.02855 + lambda22 * 199025e-9) + phi2 * (0.0998909 + phi2 * -0.0491032)) - y, deltaxDeltaLambda = 0.975534 - phi2 * (0.119161 + 3 * lambda22 * 0.0143059 + phi2 * 0.0547009), deltaxDeltaPhi = -lambdaPhi * (2 * 0.119161 + 4 * 0.0547009 * phi2 + 2 * 0.0143059 * lambda22), deltayDeltaLambda = lambdaPhi * (2 * 0.0802894 + 4 * 199025e-9 * lambda22 + 2 * -0.02855 * phi2), deltayDeltaPhi = 1.00384 + lambda22 * (0.0802894 + 199025e-9 * lambda22) + phi2 * (3 * (0.0998909 - 0.02855 * lambda22) - 5 * 0.0491032 * phi2), denominator = deltaxDeltaPhi * deltayDeltaLambda - deltayDeltaPhi * deltaxDeltaLambda, deltaLambda = (fy * deltaxDeltaPhi - fx * deltayDeltaPhi) / denominator, deltaPhi = (fx * deltayDeltaLambda - fy * deltaxDeltaLambda) / denominator;
    lambda -= deltaLambda, phi -= deltaPhi;
  } while ((abs2(deltaLambda) > epsilon3 || abs2(deltaPhi) > epsilon3) && --i > 0);
  return i && [lambda, phi];
};
function laskowski_default() {
  return projection(laskowskiRaw).scale(139.98);
}

// node_modules/d3-geo-projection/src/littrow.js
function littrowRaw(lambda, phi) {
  return [
    sin2(lambda) / cos2(phi),
    tan2(phi) * cos2(lambda)
  ];
}
littrowRaw.invert = function(x, y) {
  var x2 = x * x, y2 = y * y, y2_1 = y2 + 1, x2_y2_1 = x2 + y2_1, cosPhi = x ? sqrt1_2 * sqrt3((x2_y2_1 - sqrt3(x2_y2_1 * x2_y2_1 - 4 * x2)) / x2) : 1 / sqrt3(y2_1);
  return [
    asin2(x * cosPhi),
    sign2(y) * acos2(cosPhi)
  ];
};
function littrow_default() {
  return projection(littrowRaw).scale(144.049).clipAngle(90 - 1e-3);
}

// node_modules/d3-geo-projection/src/loximuthal.js
function loximuthalRaw(phi03) {
  var cosPhi02 = cos2(phi03), tanPhi0 = tan2(quarterPi2 + phi03 / 2);
  function forward(lambda, phi) {
    var y = phi - phi03, x = abs2(y) < epsilon3 ? lambda * cosPhi02 : abs2(x = quarterPi2 + phi / 2) < epsilon3 || abs2(abs2(x) - halfPi2) < epsilon3 ? 0 : lambda * y / log2(tan2(x) / tanPhi0);
    return [x, y];
  }
  forward.invert = function(x, y) {
    var lambda, phi = y + phi03;
    return [
      abs2(y) < epsilon3 ? x / cosPhi02 : abs2(lambda = quarterPi2 + phi / 2) < epsilon3 || abs2(abs2(lambda) - halfPi2) < epsilon3 ? 0 : x * log2(tan2(lambda) / tanPhi0) / y,
      phi
    ];
  };
  return forward;
}
function loximuthal_default() {
  return parallel1_default(loximuthalRaw).parallel(40).scale(158.837);
}

// node_modules/d3-geo-projection/src/miller.js
function millerRaw(lambda, phi) {
  return [lambda, 1.25 * log2(tan2(quarterPi2 + 0.4 * phi))];
}
millerRaw.invert = function(x, y) {
  return [x, 2.5 * atan3(exp2(0.8 * y)) - 0.625 * pi2];
};
function miller_default() {
  return projection(millerRaw).scale(108.318);
}

// node_modules/d3-geo-projection/src/mtFlatPolarParabolic.js
var sqrt6 = sqrt3(6);
var sqrt7 = sqrt3(7);
function mtFlatPolarParabolicRaw(lambda, phi) {
  var theta = asin2(7 * sin2(phi) / (3 * sqrt6));
  return [
    sqrt6 * lambda * (2 * cos2(2 * theta / 3) - 1) / sqrt7,
    9 * sin2(theta / 3) / sqrt7
  ];
}
mtFlatPolarParabolicRaw.invert = function(x, y) {
  var theta = 3 * asin2(y * sqrt7 / 9);
  return [
    x * sqrt7 / (sqrt6 * (2 * cos2(2 * theta / 3) - 1)),
    asin2(sin2(theta) * 3 * sqrt6 / 7)
  ];
};
function mtFlatPolarParabolic_default() {
  return projection(mtFlatPolarParabolicRaw).scale(164.859);
}

// node_modules/d3-geo-projection/src/mtFlatPolarQuartic.js
function mtFlatPolarQuarticRaw(lambda, phi) {
  var k2 = (1 + sqrt1_2) * sin2(phi), theta = phi;
  for (var i = 0, delta; i < 25; i++) {
    theta -= delta = (sin2(theta / 2) + sin2(theta) - k2) / (0.5 * cos2(theta / 2) + cos2(theta));
    if (abs2(delta) < epsilon3)
      break;
  }
  return [
    lambda * (1 + 2 * cos2(theta) / cos2(theta / 2)) / (3 * sqrt2),
    2 * sqrt3(3) * sin2(theta / 2) / sqrt3(2 + sqrt2)
  ];
}
mtFlatPolarQuarticRaw.invert = function(x, y) {
  var sinTheta_2 = y * sqrt3(2 + sqrt2) / (2 * sqrt3(3)), theta = 2 * asin2(sinTheta_2);
  return [
    3 * sqrt2 * x / (1 + 2 * cos2(theta) / cos2(theta / 2)),
    asin2((sinTheta_2 + sin2(theta)) / (1 + sqrt1_2))
  ];
};
function mtFlatPolarQuartic_default() {
  return projection(mtFlatPolarQuarticRaw).scale(188.209);
}

// node_modules/d3-geo-projection/src/mtFlatPolarSinusoidal.js
function mtFlatPolarSinusoidalRaw(lambda, phi) {
  var A5 = sqrt3(6 / (4 + pi2)), k2 = (1 + pi2 / 4) * sin2(phi), theta = phi / 2;
  for (var i = 0, delta; i < 25; i++) {
    theta -= delta = (theta / 2 + sin2(theta) - k2) / (0.5 + cos2(theta));
    if (abs2(delta) < epsilon3)
      break;
  }
  return [
    A5 * (0.5 + cos2(theta)) * lambda / 1.5,
    A5 * theta
  ];
}
mtFlatPolarSinusoidalRaw.invert = function(x, y) {
  var A5 = sqrt3(6 / (4 + pi2)), theta = y / A5;
  if (abs2(abs2(theta) - halfPi2) < epsilon3)
    theta = theta < 0 ? -halfPi2 : halfPi2;
  return [
    1.5 * x / (A5 * (0.5 + cos2(theta))),
    asin2((theta / 2 + sin2(theta)) / (1 + pi2 / 4))
  ];
};
function mtFlatPolarSinusoidal_default() {
  return projection(mtFlatPolarSinusoidalRaw).scale(166.518);
}

// node_modules/d3-geo-projection/src/nellHammer.js
function nellHammerRaw(lambda, phi) {
  return [
    lambda * (1 + cos2(phi)) / 2,
    2 * (phi - tan2(phi / 2))
  ];
}
nellHammerRaw.invert = function(x, y) {
  var p = y / 2;
  for (var i = 0, delta = Infinity; i < 10 && abs2(delta) > epsilon3; ++i) {
    var c = cos2(y / 2);
    y -= delta = (y - tan2(y / 2) - p) / (1 - 0.5 / (c * c));
  }
  return [
    2 * x / (1 + cos2(y)),
    y
  ];
};
function nellHammer_default() {
  return projection(nellHammerRaw).scale(152.63);
}

// node_modules/d3-geo-projection/src/patterson.js
var pattersonK1 = 1.0148;
var pattersonK2 = 0.23185;
var pattersonK3 = -0.14499;
var pattersonK4 = 0.02406;
var pattersonC1 = pattersonK1;
var pattersonC2 = 5 * pattersonK2;
var pattersonC3 = 7 * pattersonK3;
var pattersonC4 = 9 * pattersonK4;
var pattersonYmax = 1.790857183;
function pattersonRaw(lambda, phi) {
  var phi2 = phi * phi;
  return [
    lambda,
    phi * (pattersonK1 + phi2 * phi2 * (pattersonK2 + phi2 * (pattersonK3 + pattersonK4 * phi2)))
  ];
}
pattersonRaw.invert = function(x, y) {
  if (y > pattersonYmax)
    y = pattersonYmax;
  else if (y < -pattersonYmax)
    y = -pattersonYmax;
  var yc = y, delta;
  do {
    var y2 = yc * yc;
    yc -= delta = (yc * (pattersonK1 + y2 * y2 * (pattersonK2 + y2 * (pattersonK3 + pattersonK4 * y2))) - y) / (pattersonC1 + y2 * y2 * (pattersonC2 + y2 * (pattersonC3 + pattersonC4 * y2)));
  } while (abs2(delta) > epsilon3);
  return [x, yc];
};
function patterson_default() {
  return projection(pattersonRaw).scale(139.319);
}

// node_modules/d3-geo-projection/src/polyconic.js
function polyconicRaw(lambda, phi) {
  if (abs2(phi) < epsilon3)
    return [lambda, 0];
  var tanPhi = tan2(phi), k2 = lambda * sin2(phi);
  return [
    sin2(k2) / tanPhi,
    phi + (1 - cos2(k2)) / tanPhi
  ];
}
polyconicRaw.invert = function(x, y) {
  if (abs2(y) < epsilon3)
    return [x, 0];
  var k2 = x * x + y * y, phi = y * 0.5, i = 10, delta;
  do {
    var tanPhi = tan2(phi), secPhi = 1 / cos2(phi), j = k2 - 2 * y * phi + phi * phi;
    phi -= delta = (tanPhi * j + 2 * (phi - y)) / (2 + j * secPhi * secPhi + 2 * (phi - y) * tanPhi);
  } while (abs2(delta) > epsilon3 && --i > 0);
  tanPhi = tan2(phi);
  return [
    (abs2(y) < abs2(phi + 1 / tanPhi) ? asin2(x * tanPhi) : sign2(y) * sign2(x) * (acos2(abs2(x * tanPhi)) + halfPi2)) / sin2(phi),
    phi
  ];
};
function polyconic_default() {
  return projection(polyconicRaw).scale(103.74);
}

// node_modules/d3-geo-projection/src/polyhedral/matrix.js
function matrix_default(a, b) {
  var u = subtract(a[1], a[0]), v = subtract(b[1], b[0]), phi = angle2(u, v), s = length(u) / length(v);
  return multiply([
    1,
    0,
    a[0][0],
    0,
    1,
    a[0][1]
  ], multiply([
    s,
    0,
    0,
    0,
    s,
    0
  ], multiply([
    cos2(phi),
    sin2(phi),
    0,
    -sin2(phi),
    cos2(phi),
    0
  ], [
    1,
    0,
    -b[0][0],
    0,
    1,
    -b[0][1]
  ])));
}
function inverse(m) {
  var k2 = 1 / (m[0] * m[4] - m[1] * m[3]);
  return [
    k2 * m[4],
    -k2 * m[1],
    k2 * (m[1] * m[5] - m[2] * m[4]),
    -k2 * m[3],
    k2 * m[0],
    k2 * (m[2] * m[3] - m[0] * m[5])
  ];
}
function multiply(a, b) {
  return [
    a[0] * b[0] + a[1] * b[3],
    a[0] * b[1] + a[1] * b[4],
    a[0] * b[2] + a[1] * b[5] + a[2],
    a[3] * b[0] + a[4] * b[3],
    a[3] * b[1] + a[4] * b[4],
    a[3] * b[2] + a[4] * b[5] + a[5]
  ];
}
function subtract(a, b) {
  return [a[0] - b[0], a[1] - b[1]];
}
function length(v) {
  return sqrt3(v[0] * v[0] + v[1] * v[1]);
}
function angle2(a, b) {
  return atan22(a[0] * b[1] - a[1] * b[0], a[0] * b[0] + a[1] * b[1]);
}

// node_modules/d3-geo-projection/src/polyhedral/index.js
function polyhedral_default(root, face) {
  recurse(root, { transform: null });
  function recurse(node, parent) {
    node.edges = faceEdges(node.face);
    if (parent.face) {
      var shared = node.shared = sharedEdge(node.face, parent.face), m = matrix_default(shared.map(parent.project), shared.map(node.project));
      node.transform = parent.transform ? multiply(parent.transform, m) : m;
      var edges = parent.edges;
      for (var i = 0, n = edges.length; i < n; ++i) {
        if (pointEqual2(shared[0], edges[i][1]) && pointEqual2(shared[1], edges[i][0]))
          edges[i] = node;
        if (pointEqual2(shared[0], edges[i][0]) && pointEqual2(shared[1], edges[i][1]))
          edges[i] = node;
      }
      edges = node.edges;
      for (i = 0, n = edges.length; i < n; ++i) {
        if (pointEqual2(shared[0], edges[i][0]) && pointEqual2(shared[1], edges[i][1]))
          edges[i] = parent;
        if (pointEqual2(shared[0], edges[i][1]) && pointEqual2(shared[1], edges[i][0]))
          edges[i] = parent;
      }
    } else {
      node.transform = parent.transform;
    }
    if (node.children) {
      node.children.forEach(function(child) {
        recurse(child, node);
      });
    }
    return node;
  }
  function forward(lambda, phi) {
    var node = face(lambda, phi), point = node.project([lambda * degrees2, phi * degrees2]), t;
    if (t = node.transform) {
      return [
        t[0] * point[0] + t[1] * point[1] + t[2],
        -(t[3] * point[0] + t[4] * point[1] + t[5])
      ];
    }
    point[1] = -point[1];
    return point;
  }
  if (hasInverse(root))
    forward.invert = function(x, y) {
      var coordinates = faceInvert(root, [x, -y]);
      return coordinates && (coordinates[0] *= radians2, coordinates[1] *= radians2, coordinates);
    };
  function faceInvert(node, coordinates) {
    var invert = node.project.invert, t = node.transform, point = coordinates;
    if (t) {
      t = inverse(t);
      point = [
        t[0] * point[0] + t[1] * point[1] + t[2],
        t[3] * point[0] + t[4] * point[1] + t[5]
      ];
    }
    if (invert && node === faceDegrees(p = invert(point)))
      return p;
    var p, children = node.children;
    for (var i = 0, n = children && children.length; i < n; ++i) {
      if (p = faceInvert(children[i], coordinates))
        return p;
    }
  }
  function faceDegrees(coordinates) {
    return face(coordinates[0] * radians2, coordinates[1] * radians2);
  }
  var proj = projection(forward), stream_ = proj.stream;
  proj.stream = function(stream) {
    var rotate = proj.rotate(), rotateStream = stream_(stream), sphereStream = (proj.rotate([0, 0]), stream_(stream));
    proj.rotate(rotate);
    rotateStream.sphere = function() {
      sphereStream.polygonStart();
      sphereStream.lineStart();
      outline(sphereStream, root);
      sphereStream.lineEnd();
      sphereStream.polygonEnd();
    };
    return rotateStream;
  };
  return proj.angle(-30);
}
function outline(stream, node, parent) {
  var point, edges = node.edges, n = edges.length, edge, multiPoint = { type: "MultiPoint", coordinates: node.face }, notPoles = node.face.filter(function(d) {
    return abs2(d[1]) !== 90;
  }), b = bounds_default({ type: "MultiPoint", coordinates: notPoles }), inside = false, j = -1, dx = b[1][0] - b[0][0];
  var c = dx === 180 || dx === 360 ? [(b[0][0] + b[1][0]) / 2, (b[0][1] + b[1][1]) / 2] : centroid_default(multiPoint);
  if (parent)
    while (++j < n) {
      if (edges[j] === parent)
        break;
    }
  ++j;
  for (var i = 0; i < n; ++i) {
    edge = edges[(i + j) % n];
    if (Array.isArray(edge)) {
      if (!inside) {
        stream.point((point = interpolate_default(edge[0], c)(epsilon3))[0], point[1]);
        inside = true;
      }
      stream.point((point = interpolate_default(edge[1], c)(epsilon3))[0], point[1]);
    } else {
      inside = false;
      if (edge !== parent)
        outline(stream, edge, node);
    }
  }
}
function pointEqual2(a, b) {
  return a && b && a[0] === b[0] && a[1] === b[1];
}
function sharedEdge(a, b) {
  var x, y, n = a.length, found = null;
  for (var i = 0; i < n; ++i) {
    x = a[i];
    for (var j = b.length; --j >= 0; ) {
      y = b[j];
      if (x[0] === y[0] && x[1] === y[1]) {
        if (found)
          return [found, x];
        found = x;
      }
    }
  }
}
function faceEdges(face) {
  var n = face.length, edges = [];
  for (var a = face[n - 1], i = 0; i < n; ++i)
    edges.push([a, a = face[i]]);
  return edges;
}
function hasInverse(node) {
  return node.project.invert || node.children && node.children.some(hasInverse);
}

// node_modules/d3-geo-projection/src/polyhedral/octahedron.js
var octahedron = [
  [0, 90],
  [-90, 0],
  [0, 0],
  [90, 0],
  [180, 0],
  [0, -90]
];
var octahedron_default = [
  [0, 2, 1],
  [0, 3, 2],
  [5, 1, 2],
  [5, 2, 3],
  [0, 1, 4],
  [0, 4, 3],
  [5, 4, 1],
  [5, 3, 4]
].map(function(face) {
  return face.map(function(i) {
    return octahedron[i];
  });
});

// node_modules/d3-geo-projection/src/polyhedral/butterfly.js
function butterfly_default(faceProjection) {
  faceProjection = faceProjection || function(face) {
    var c = centroid_default({ type: "MultiPoint", coordinates: face });
    return gnomonic_default().scale(1).translate([0, 0]).rotate([-c[0], -c[1]]);
  };
  var faces = octahedron_default.map(function(face) {
    return { face, project: faceProjection(face) };
  });
  [-1, 0, 0, 1, 0, 1, 4, 5].forEach(function(d, i) {
    var node = faces[d];
    node && (node.children || (node.children = [])).push(faces[i]);
  });
  return polyhedral_default(faces[0], function(lambda, phi) {
    return faces[lambda < -pi2 / 2 ? phi < 0 ? 6 : 4 : lambda < 0 ? phi < 0 ? 2 : 0 : lambda < pi2 / 2 ? phi < 0 ? 3 : 1 : phi < 0 ? 7 : 5];
  }).angle(-30).scale(101.858).center([0, 45]);
}

// node_modules/d3-geo-projection/src/polyhedral/collignon.js
var kx = 2 / sqrt3(3);
function collignonK(a, b) {
  var p = collignonRaw(a, b);
  return [p[0] * kx, p[1]];
}
collignonK.invert = function(x, y) {
  return collignonRaw.invert(x / kx, y);
};
function collignon_default2(faceProjection) {
  faceProjection = faceProjection || function(face) {
    var c = centroid_default({ type: "MultiPoint", coordinates: face });
    return projection(collignonK).translate([0, 0]).scale(1).rotate(c[1] > 0 ? [-c[0], 0] : [180 - c[0], 180]);
  };
  var faces = octahedron_default.map(function(face) {
    return { face, project: faceProjection(face) };
  });
  [-1, 0, 0, 1, 0, 1, 4, 5].forEach(function(d, i) {
    var node = faces[d];
    node && (node.children || (node.children = [])).push(faces[i]);
  });
  return polyhedral_default(faces[0], function(lambda, phi) {
    return faces[lambda < -pi2 / 2 ? phi < 0 ? 6 : 4 : lambda < 0 ? phi < 0 ? 2 : 0 : lambda < pi2 / 2 ? phi < 0 ? 3 : 1 : phi < 0 ? 7 : 5];
  }).angle(-30).scale(121.906).center([0, 48.5904]);
}

// node_modules/d3-geo-projection/src/polyhedral/waterman.js
function waterman_default(faceProjection) {
  faceProjection = faceProjection || function(face2) {
    var c = face2.length === 6 ? centroid_default({ type: "MultiPoint", coordinates: face2 }) : face2[0];
    return gnomonic_default().scale(1).translate([0, 0]).rotate([-c[0], -c[1]]);
  };
  var w5 = octahedron_default.map(function(face2) {
    var xyz = face2.map(cartesian2), n = xyz.length, a = xyz[n - 1], b, hexagon = [];
    for (var i = 0; i < n; ++i) {
      b = xyz[i];
      hexagon.push(spherical2([
        a[0] * 0.9486832980505138 + b[0] * 0.31622776601683794,
        a[1] * 0.9486832980505138 + b[1] * 0.31622776601683794,
        a[2] * 0.9486832980505138 + b[2] * 0.31622776601683794
      ]), spherical2([
        b[0] * 0.9486832980505138 + a[0] * 0.31622776601683794,
        b[1] * 0.9486832980505138 + a[1] * 0.31622776601683794,
        b[2] * 0.9486832980505138 + a[2] * 0.31622776601683794
      ]));
      a = b;
    }
    return hexagon;
  });
  var cornerNormals = [];
  var parents = [-1, 0, 0, 1, 0, 1, 4, 5];
  w5.forEach(function(hexagon, j) {
    var face2 = octahedron_default[j], n = face2.length, normals = cornerNormals[j] = [];
    for (var i = 0; i < n; ++i) {
      w5.push([
        face2[i],
        hexagon[(i * 2 + 2) % (2 * n)],
        hexagon[(i * 2 + 1) % (2 * n)]
      ]);
      parents.push(j);
      normals.push(cross(
        cartesian2(hexagon[(i * 2 + 2) % (2 * n)]),
        cartesian2(hexagon[(i * 2 + 1) % (2 * n)])
      ));
    }
  });
  var faces = w5.map(function(face2) {
    return {
      project: faceProjection(face2),
      face: face2
    };
  });
  parents.forEach(function(d, i) {
    var parent = faces[d];
    parent && (parent.children || (parent.children = [])).push(faces[i]);
  });
  function face(lambda, phi) {
    var cosphi = cos2(phi), p = [cosphi * cos2(lambda), cosphi * sin2(lambda), sin2(phi)];
    var hexagon = lambda < -pi2 / 2 ? phi < 0 ? 6 : 4 : lambda < 0 ? phi < 0 ? 2 : 0 : lambda < pi2 / 2 ? phi < 0 ? 3 : 1 : phi < 0 ? 7 : 5;
    var n = cornerNormals[hexagon];
    return faces[dot(n[0], p) < 0 ? 8 + 3 * hexagon : dot(n[1], p) < 0 ? 8 + 3 * hexagon + 1 : dot(n[2], p) < 0 ? 8 + 3 * hexagon + 2 : hexagon];
  }
  return polyhedral_default(faces[0], face).angle(-30).scale(110.625).center([0, 45]);
}
function dot(a, b) {
  for (var i = 0, n = a.length, s = 0; i < n; ++i)
    s += a[i] * b[i];
  return s;
}
function cross(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
  ];
}
function spherical2(cartesian3) {
  return [
    atan22(cartesian3[1], cartesian3[0]) * degrees2,
    asin2(max(-1, min(1, cartesian3[2]))) * degrees2
  ];
}
function cartesian2(coordinates) {
  var lambda = coordinates[0] * radians2, phi = coordinates[1] * radians2, cosphi = cos2(phi);
  return [
    cosphi * cos2(lambda),
    cosphi * sin2(lambda),
    sin2(phi)
  ];
}

// node_modules/d3-geo-projection/src/noop.js
var noop_default = () => {
};

// node_modules/d3-geo-projection/src/project/clockwise.js
function clockwise_default(ring) {
  if ((n = ring.length) < 4)
    return false;
  var i = 0, n, area = ring[n - 1][1] * ring[0][0] - ring[n - 1][0] * ring[0][1];
  while (++i < n)
    area += ring[i - 1][1] * ring[i][0] - ring[i - 1][0] * ring[i][1];
  return area <= 0;
}

// node_modules/d3-geo-projection/src/project/contains.js
function contains_default(ring, point) {
  var x = point[0], y = point[1], contains = false;
  for (var i = 0, n = ring.length, j = n - 1; i < n; j = i++) {
    var pi3 = ring[i], xi = pi3[0], yi = pi3[1], pj = ring[j], xj = pj[0], yj = pj[1];
    if (yi > y ^ yj > y && x < (xj - xi) * (y - yi) / (yj - yi) + xi)
      contains = !contains;
  }
  return contains;
}

// node_modules/d3-geo-projection/src/project/index.js
function project_default(object, projection2) {
  var stream = projection2.stream, project;
  if (!stream)
    throw new Error("invalid projection");
  switch (object && object.type) {
    case "Feature":
      project = projectFeature;
      break;
    case "FeatureCollection":
      project = projectFeatureCollection;
      break;
    default:
      project = projectGeometry;
      break;
  }
  return project(object, stream);
}
function projectFeatureCollection(o, stream) {
  return {
    type: "FeatureCollection",
    features: o.features.map(function(f) {
      return projectFeature(f, stream);
    })
  };
}
function projectFeature(o, stream) {
  return {
    type: "Feature",
    id: o.id,
    properties: o.properties,
    geometry: projectGeometry(o.geometry, stream)
  };
}
function projectGeometryCollection(o, stream) {
  return {
    type: "GeometryCollection",
    geometries: o.geometries.map(function(o2) {
      return projectGeometry(o2, stream);
    })
  };
}
function projectGeometry(o, stream) {
  if (!o)
    return null;
  if (o.type === "GeometryCollection")
    return projectGeometryCollection(o, stream);
  var sink;
  switch (o.type) {
    case "Point":
      sink = sinkPoint;
      break;
    case "MultiPoint":
      sink = sinkPoint;
      break;
    case "LineString":
      sink = sinkLine;
      break;
    case "MultiLineString":
      sink = sinkLine;
      break;
    case "Polygon":
      sink = sinkPolygon;
      break;
    case "MultiPolygon":
      sink = sinkPolygon;
      break;
    case "Sphere":
      sink = sinkPolygon;
      break;
    default:
      return null;
  }
  stream_default(o, stream(sink));
  return sink.result();
}
var points = [];
var lines = [];
var sinkPoint = {
  point: function(x, y) {
    points.push([x, y]);
  },
  result: function() {
    var result = !points.length ? null : points.length < 2 ? { type: "Point", coordinates: points[0] } : { type: "MultiPoint", coordinates: points };
    points = [];
    return result;
  }
};
var sinkLine = {
  lineStart: noop_default,
  point: function(x, y) {
    points.push([x, y]);
  },
  lineEnd: function() {
    if (points.length)
      lines.push(points), points = [];
  },
  result: function() {
    var result = !lines.length ? null : lines.length < 2 ? { type: "LineString", coordinates: lines[0] } : { type: "MultiLineString", coordinates: lines };
    lines = [];
    return result;
  }
};
var sinkPolygon = {
  polygonStart: noop_default,
  lineStart: noop_default,
  point: function(x, y) {
    points.push([x, y]);
  },
  lineEnd: function() {
    var n = points.length;
    if (n) {
      do
        points.push(points[0].slice());
      while (++n < 4);
      lines.push(points), points = [];
    }
  },
  polygonEnd: noop_default,
  result: function() {
    if (!lines.length)
      return null;
    var polygons = [], holes = [];
    lines.forEach(function(ring) {
      if (clockwise_default(ring))
        polygons.push([ring]);
      else
        holes.push(ring);
    });
    holes.forEach(function(hole) {
      var point = hole[0];
      polygons.some(function(polygon) {
        if (contains_default(polygon[0], point)) {
          polygon.push(hole);
          return true;
        }
      }) || polygons.push([hole]);
    });
    lines = [];
    return !polygons.length ? null : polygons.length > 1 ? { type: "MultiPolygon", coordinates: polygons } : { type: "Polygon", coordinates: polygons[0] };
  }
};

// node_modules/d3-geo-projection/src/quincuncial/index.js
function quincuncial_default(project) {
  var dx = project(halfPi2, 0)[0] - project(-halfPi2, 0)[0];
  function projectQuincuncial(lambda, phi) {
    var t = abs2(lambda) < halfPi2, p = project(t ? lambda : lambda > 0 ? lambda - pi2 : lambda + pi2, phi), x = (p[0] - p[1]) * sqrt1_2, y = (p[0] + p[1]) * sqrt1_2;
    if (t)
      return [x, y];
    var d = dx * sqrt1_2, s = x > 0 ^ y > 0 ? -1 : 1;
    return [s * x - sign2(y) * d, s * y - sign2(x) * d];
  }
  if (project.invert)
    projectQuincuncial.invert = function(x03, y03) {
      var x = (x03 + y03) * sqrt1_2, y = (y03 - x03) * sqrt1_2, t = abs2(x) < 0.5 * dx && abs2(y) < 0.5 * dx;
      if (!t) {
        var d = dx * sqrt1_2, s = x > 0 ^ y > 0 ? -1 : 1, x12 = -s * x03 + (y > 0 ? 1 : -1) * d, y12 = -s * y03 + (x > 0 ? 1 : -1) * d;
        x = (-x12 - y12) * sqrt1_2;
        y = (x12 - y12) * sqrt1_2;
      }
      var p = project.invert(x, y);
      if (!t)
        p[0] += x > 0 ? pi2 : -pi2;
      return p;
    };
  return projection(projectQuincuncial).rotate([-90, -90, 45]).clipAngle(180 - 1e-3);
}

// node_modules/d3-geo-projection/src/quincuncial/gringorten.js
function gringorten_default2() {
  return quincuncial_default(gringortenRaw).scale(176.423);
}

// node_modules/d3-geo-projection/src/quincuncial/peirce.js
function peirce_default() {
  return quincuncial_default(guyouRaw).scale(111.48);
}

// node_modules/d3-geo-projection/src/rectangularPolyconic.js
function rectangularPolyconicRaw(phi03) {
  var sinPhi02 = sin2(phi03);
  function forward(lambda, phi) {
    var A5 = sinPhi02 ? tan2(lambda * sinPhi02 / 2) / sinPhi02 : lambda / 2;
    if (!phi)
      return [2 * A5, -phi03];
    var E = 2 * atan3(A5 * sin2(phi)), cotPhi = 1 / tan2(phi);
    return [
      sin2(E) * cotPhi,
      phi + (1 - cos2(E)) * cotPhi - phi03
    ];
  }
  forward.invert = function(x, y) {
    if (abs2(y += phi03) < epsilon3)
      return [sinPhi02 ? 2 * atan3(sinPhi02 * x / 2) / sinPhi02 : x, 0];
    var k2 = x * x + y * y, phi = 0, i = 10, delta;
    do {
      var tanPhi = tan2(phi), secPhi = 1 / cos2(phi), j = k2 - 2 * y * phi + phi * phi;
      phi -= delta = (tanPhi * j + 2 * (phi - y)) / (2 + j * secPhi * secPhi + 2 * (phi - y) * tanPhi);
    } while (abs2(delta) > epsilon3 && --i > 0);
    var E = x * (tanPhi = tan2(phi)), A5 = tan2(abs2(y) < abs2(phi + 1 / tanPhi) ? asin2(E) * 0.5 : acos2(E) * 0.5 + pi2 / 4) / sin2(phi);
    return [
      sinPhi02 ? 2 * atan3(sinPhi02 * A5) / sinPhi02 : 2 * A5,
      phi
    ];
  };
  return forward;
}
function rectangularPolyconic_default() {
  return parallel1_default(rectangularPolyconicRaw).scale(131.215);
}

// node_modules/d3-geo-projection/src/robinson.js
var K2 = [
  [0.9986, -0.062],
  [1, 0],
  [0.9986, 0.062],
  [0.9954, 0.124],
  [0.99, 0.186],
  [0.9822, 0.248],
  [0.973, 0.31],
  [0.96, 0.372],
  [0.9427, 0.434],
  [0.9216, 0.4958],
  [0.8962, 0.5571],
  [0.8679, 0.6176],
  [0.835, 0.6769],
  [0.7986, 0.7346],
  [0.7597, 0.7903],
  [0.7186, 0.8435],
  [0.6732, 0.8936],
  [0.6213, 0.9394],
  [0.5722, 0.9761],
  [0.5322, 1]
];
K2.forEach(function(d) {
  d[1] *= 1.593415793900743;
});
function robinsonRaw(lambda, phi) {
  var i = min(18, abs2(phi) * 36 / pi2), i0 = floor(i), di = i - i0, ax = (k2 = K2[i0])[0], ay = k2[1], bx = (k2 = K2[++i0])[0], by = k2[1], cx = (k2 = K2[min(19, ++i0)])[0], cy = k2[1], k2;
  return [
    lambda * (bx + di * (cx - ax) / 2 + di * di * (cx - 2 * bx + ax) / 2),
    sign2(phi) * (by + di * (cy - ay) / 2 + di * di * (cy - 2 * by + ay) / 2)
  ];
}
robinsonRaw.invert = function(x, y) {
  var phi = y * 90, i = min(18, abs2(phi / 5)), i0 = max(0, floor(i));
  do {
    var ay = K2[i0][1], by = K2[i0 + 1][1], cy = K2[min(19, i0 + 2)][1], u = cy - ay, v = cy - 2 * by + ay, t = 2 * (abs2(y) - by) / u, c = v / u, di = t * (1 - c * t * (1 - 2 * c * t));
    if (di >= 0 || i0 === 1) {
      phi = (y >= 0 ? 5 : -5) * (di + i);
      var j = 50, delta;
      do {
        i = min(18, abs2(phi) / 5);
        i0 = floor(i);
        di = i - i0;
        ay = K2[i0][1];
        by = K2[i0 + 1][1];
        cy = K2[min(19, i0 + 2)][1];
        phi -= (delta = sign2(y) * (by + di * (cy - ay) / 2 + di * di * (cy - 2 * by + ay) / 2) - y) * degrees2;
      } while (abs2(delta) > epsilon22 && --j > 0);
      break;
    }
  } while (--i0 >= 0);
  var ax = K2[i0][0], bx = K2[i0 + 1][0], cx = K2[min(19, i0 + 2)][0];
  return [
    x / (bx + di * (cx - ax) / 2 + di * di * (cx - 2 * bx + ax) / 2),
    phi * radians2
  ];
};
function robinson_default() {
  return projection(robinsonRaw).scale(152.63);
}

// node_modules/d3-geo-projection/src/satellite.js
function satelliteVerticalRaw(P) {
  function forward(lambda, phi) {
    var cosPhi = cos2(phi), k2 = (P - 1) / (P - cosPhi * cos2(lambda));
    return [
      k2 * cosPhi * sin2(lambda),
      k2 * sin2(phi)
    ];
  }
  forward.invert = function(x, y) {
    var rho2 = x * x + y * y, rho = sqrt3(rho2), sinc = (P - sqrt3(1 - rho2 * (P + 1) / (P - 1))) / ((P - 1) / rho + rho / (P - 1));
    return [
      atan22(x * sinc, rho * sqrt3(1 - sinc * sinc)),
      rho ? asin2(y * sinc / rho) : 0
    ];
  };
  return forward;
}
function satelliteRaw(P, omega) {
  var vertical = satelliteVerticalRaw(P);
  if (!omega)
    return vertical;
  var cosOmega = cos2(omega), sinOmega = sin2(omega);
  function forward(lambda, phi) {
    var coordinates = vertical(lambda, phi), y = coordinates[1], A5 = y * sinOmega / (P - 1) + cosOmega;
    return [
      coordinates[0] * cosOmega / A5,
      y / A5
    ];
  }
  forward.invert = function(x, y) {
    var k2 = (P - 1) / (P - 1 - y * sinOmega);
    return vertical.invert(k2 * x, k2 * y * cosOmega);
  };
  return forward;
}
function satellite_default() {
  var distance = 2, omega = 0, m = projectionMutator(satelliteRaw), p = m(distance, omega);
  p.distance = function(_) {
    if (!arguments.length)
      return distance;
    return m(distance = +_, omega);
  };
  p.tilt = function(_) {
    if (!arguments.length)
      return omega * degrees2;
    return m(distance, omega = _ * radians2);
  };
  return p.scale(432.147).clipAngle(acos2(1 / distance) * degrees2 - 1e-6);
}

// node_modules/d3-geo-projection/src/times.js
function timesRaw(lambda, phi) {
  var t = tan2(phi / 2), s = sin2(quarterPi2 * t);
  return [
    lambda * (0.74482 - 0.34588 * s * s),
    1.70711 * t
  ];
}
timesRaw.invert = function(x, y) {
  var t = y / 1.70711, s = sin2(quarterPi2 * t);
  return [
    x / (0.74482 - 0.34588 * s * s),
    2 * atan3(t)
  ];
};
function times_default() {
  return projection(timesRaw).scale(146.153);
}

// node_modules/d3-geo-projection/src/vanDerGrinten.js
function vanDerGrintenRaw(lambda, phi) {
  if (abs2(phi) < epsilon3)
    return [lambda, 0];
  var sinTheta = abs2(phi / halfPi2), theta = asin2(sinTheta);
  if (abs2(lambda) < epsilon3 || abs2(abs2(phi) - halfPi2) < epsilon3)
    return [0, sign2(phi) * pi2 * tan2(theta / 2)];
  var cosTheta = cos2(theta), A5 = abs2(pi2 / lambda - lambda / pi2) / 2, A22 = A5 * A5, G = cosTheta / (sinTheta + cosTheta - 1), P = G * (2 / sinTheta - 1), P2 = P * P, P2_A2 = P2 + A22, G_P2 = G - P2, Q = A22 + G;
  return [
    sign2(lambda) * pi2 * (A5 * G_P2 + sqrt3(A22 * G_P2 * G_P2 - P2_A2 * (G * G - P2))) / P2_A2,
    sign2(phi) * pi2 * (P * Q - A5 * sqrt3((A22 + 1) * P2_A2 - Q * Q)) / P2_A2
  ];
}
vanDerGrintenRaw.invert = function(x, y) {
  if (abs2(y) < epsilon3)
    return [x, 0];
  if (abs2(x) < epsilon3)
    return [0, halfPi2 * sin2(2 * atan3(y / pi2))];
  var x2 = (x /= pi2) * x, y2 = (y /= pi2) * y, x2_y2 = x2 + y2, z = x2_y2 * x2_y2, c1 = -abs2(y) * (1 + x2_y2), c2 = c1 - 2 * y2 + x2, c3 = -2 * c1 + 1 + 2 * y2 + z, d = y2 / c3 + (2 * c2 * c2 * c2 / (c3 * c3 * c3) - 9 * c1 * c2 / (c3 * c3)) / 27, a1 = (c1 - c2 * c2 / (3 * c3)) / c3, m1 = 2 * sqrt3(-a1 / 3), theta1 = acos2(3 * d / (a1 * m1)) / 3;
  return [
    pi2 * (x2_y2 - 1 + sqrt3(1 + 2 * (x2 - y2) + z)) / (2 * x),
    sign2(y) * pi2 * (-m1 * cos2(theta1 + pi2 / 3) - c2 / (3 * c3))
  ];
};
function vanDerGrinten_default() {
  return projection(vanDerGrintenRaw).scale(79.4183);
}

// node_modules/d3-geo-projection/src/vanDerGrinten2.js
function vanDerGrinten2Raw(lambda, phi) {
  if (abs2(phi) < epsilon3)
    return [lambda, 0];
  var sinTheta = abs2(phi / halfPi2), theta = asin2(sinTheta);
  if (abs2(lambda) < epsilon3 || abs2(abs2(phi) - halfPi2) < epsilon3)
    return [0, sign2(phi) * pi2 * tan2(theta / 2)];
  var cosTheta = cos2(theta), A5 = abs2(pi2 / lambda - lambda / pi2) / 2, A22 = A5 * A5, x12 = cosTheta * (sqrt3(1 + A22) - A5 * cosTheta) / (1 + A22 * sinTheta * sinTheta);
  return [
    sign2(lambda) * pi2 * x12,
    sign2(phi) * pi2 * sqrt3(1 - x12 * (2 * A5 + x12))
  ];
}
vanDerGrinten2Raw.invert = function(x, y) {
  if (!x)
    return [0, halfPi2 * sin2(2 * atan3(y / pi2))];
  var x12 = abs2(x / pi2), A5 = (1 - x12 * x12 - (y /= pi2) * y) / (2 * x12), A22 = A5 * A5, B2 = sqrt3(A22 + 1);
  return [
    sign2(x) * pi2 * (B2 - A5),
    sign2(y) * halfPi2 * sin2(2 * atan22(sqrt3((1 - 2 * A5 * x12) * (A5 + B2) - x12), sqrt3(B2 + A5 + x12)))
  ];
};
function vanDerGrinten2_default() {
  return projection(vanDerGrinten2Raw).scale(79.4183);
}

// node_modules/d3-geo-projection/src/vanDerGrinten3.js
function vanDerGrinten3Raw(lambda, phi) {
  if (abs2(phi) < epsilon3)
    return [lambda, 0];
  var sinTheta = phi / halfPi2, theta = asin2(sinTheta);
  if (abs2(lambda) < epsilon3 || abs2(abs2(phi) - halfPi2) < epsilon3)
    return [0, pi2 * tan2(theta / 2)];
  var A5 = (pi2 / lambda - lambda / pi2) / 2, y12 = sinTheta / (1 + cos2(theta));
  return [
    pi2 * (sign2(lambda) * sqrt3(A5 * A5 + 1 - y12 * y12) - A5),
    pi2 * y12
  ];
}
vanDerGrinten3Raw.invert = function(x, y) {
  if (!y)
    return [x, 0];
  var y12 = y / pi2, A5 = (pi2 * pi2 * (1 - y12 * y12) - x * x) / (2 * pi2 * x);
  return [
    x ? pi2 * (sign2(x) * sqrt3(A5 * A5 + 1) - A5) : 0,
    halfPi2 * sin2(2 * atan3(y12))
  ];
};
function vanDerGrinten3_default() {
  return projection(vanDerGrinten3Raw).scale(79.4183);
}

// node_modules/d3-geo-projection/src/vanDerGrinten4.js
function vanDerGrinten4Raw(lambda, phi) {
  if (!phi)
    return [lambda, 0];
  var phi03 = abs2(phi);
  if (!lambda || phi03 === halfPi2)
    return [0, phi];
  var B2 = phi03 / halfPi2, B22 = B2 * B2, C = (8 * B2 - B22 * (B22 + 2) - 5) / (2 * B22 * (B2 - 1)), C2 = C * C, BC = B2 * C, B_C2 = B22 + C2 + 2 * BC, B_3C = B2 + 3 * C, lambda03 = lambda / halfPi2, lambda12 = lambda03 + 1 / lambda03, D = sign2(abs2(lambda) - halfPi2) * sqrt3(lambda12 * lambda12 - 4), D2 = D * D, F = B_C2 * (B22 + C2 * D2 - 1) + (1 - B22) * (B22 * (B_3C * B_3C + 4 * C2) + 12 * BC * C2 + 4 * C2 * C2), x12 = (D * (B_C2 + C2 - 1) + 2 * sqrt3(F)) / (4 * B_C2 + D2);
  return [
    sign2(lambda) * halfPi2 * x12,
    sign2(phi) * halfPi2 * sqrt3(1 + D * abs2(x12) - x12 * x12)
  ];
}
vanDerGrinten4Raw.invert = function(x, y) {
  var delta;
  if (!x || !y)
    return [x, y];
  var sy = sign2(y);
  y = abs2(y) / pi2;
  var x12 = sign2(x) * x / halfPi2, D = (x12 * x12 - 1 + 4 * y * y) / abs2(x12), D2 = D * D, B2 = y * (2 - (y > 0.5 ? min(y, abs2(x)) : 0)), r = x * x + y * y, i = 50;
  do {
    var B22 = B2 * B2, C = (8 * B2 - B22 * (B22 + 2) - 5) / (2 * B22 * (B2 - 1)), C_ = (3 * B2 - B22 * B2 - 10) / (2 * B22 * B2), C2 = C * C, BC = B2 * C, B_C = B2 + C, B_C2 = B_C * B_C, B_3C = B2 + 3 * C, F = B_C2 * (B22 + C2 * D2 - 1) + (1 - B22) * (B22 * (B_3C * B_3C + 4 * C2) + C2 * (12 * BC + 4 * C2)), F_ = -2 * B_C * (4 * BC * C2 + (1 - 4 * B22 + 3 * B22 * B22) * (1 + C_) + C2 * (-6 + 14 * B22 - D2 + (-8 + 8 * B22 - 2 * D2) * C_) + BC * (-8 + 12 * B22 + (-10 + 10 * B22 - D2) * C_)), sqrtF = sqrt3(F), f = D * (B_C2 + C2 - 1) + 2 * sqrtF - x12 * (4 * B_C2 + D2), f_ = D * (2 * C * C_ + 2 * B_C * (1 + C_)) + F_ / sqrtF - 8 * B_C * (D * (-1 + C2 + B_C2) + 2 * sqrtF) * (1 + C_) / (D2 + 4 * B_C2);
    B2 -= delta = f / f_;
  } while (delta * r * r > epsilon3 && --i > 0);
  return [
    sign2(x) * (sqrt3(D * D + 4) + D) * pi2 / 4,
    sy * halfPi2 * B2
  ];
};
function vanDerGrinten4_default() {
  return projection(vanDerGrinten4Raw).scale(127.16);
}

// node_modules/d3-geo-projection/src/wagner4.js
var A = 4 * pi2 + 3 * sqrt3(3);
var B = 2 * sqrt3(2 * pi2 * sqrt3(3) / A);
var wagner4Raw = mollweideBromleyRaw(B * sqrt3(3) / pi2, B, A / 6);
function wagner4_default() {
  return projection(wagner4Raw).scale(176.84);
}

// node_modules/d3-geo-projection/src/wagner6.js
function wagner6Raw(lambda, phi) {
  return [lambda * sqrt3(1 - 3 * phi * phi / (pi2 * pi2)), phi];
}
wagner6Raw.invert = function(x, y) {
  return [x / sqrt3(1 - 3 * y * y / (pi2 * pi2)), y];
};
function wagner6_default() {
  return projection(wagner6Raw).scale(152.63);
}

// node_modules/d3-geo-projection/src/wiechel.js
function wiechelRaw(lambda, phi) {
  var cosPhi = cos2(phi), sinPhi = cos2(lambda) * cosPhi, sin1_Phi = 1 - sinPhi, cosLambda = cos2(lambda = atan22(sin2(lambda) * cosPhi, -sin2(phi))), sinLambda = sin2(lambda);
  cosPhi = sqrt3(1 - sinPhi * sinPhi);
  return [
    sinLambda * cosPhi - cosLambda * sin1_Phi,
    -cosLambda * cosPhi - sinLambda * sin1_Phi
  ];
}
wiechelRaw.invert = function(x, y) {
  var w2 = (x * x + y * y) / -2, k2 = sqrt3(-w2 * (2 + w2)), b = y * w2 + x * k2, a = x * w2 - y * k2, D = sqrt3(a * a + b * b);
  return [
    atan22(k2 * b, D * (1 + w2)),
    D ? -asin2(k2 * a / D) : 0
  ];
};
function wiechel_default() {
  return projection(wiechelRaw).rotate([0, -90, 45]).scale(124.75).clipAngle(180 - 1e-3);
}

// node_modules/d3-geo-projection/src/winkel3.js
function winkel3Raw(lambda, phi) {
  var coordinates = aitoffRaw(lambda, phi);
  return [
    (coordinates[0] + lambda / halfPi2) / 2,
    (coordinates[1] + phi) / 2
  ];
}
winkel3Raw.invert = function(x, y) {
  var lambda = x, phi = y, i = 25;
  do {
    var cosphi = cos2(phi), sinphi = sin2(phi), sin_2phi = sin2(2 * phi), sin2phi = sinphi * sinphi, cos2phi = cosphi * cosphi, sinlambda = sin2(lambda), coslambda_2 = cos2(lambda / 2), sinlambda_2 = sin2(lambda / 2), sin2lambda_2 = sinlambda_2 * sinlambda_2, C = 1 - cos2phi * coslambda_2 * coslambda_2, E = C ? acos2(cosphi * coslambda_2) * sqrt3(F = 1 / C) : F = 0, F, fx = 0.5 * (2 * E * cosphi * sinlambda_2 + lambda / halfPi2) - x, fy = 0.5 * (E * sinphi + phi) - y, dxdlambda = 0.5 * F * (cos2phi * sin2lambda_2 + E * cosphi * coslambda_2 * sin2phi) + 0.5 / halfPi2, dxdphi = F * (sinlambda * sin_2phi / 4 - E * sinphi * sinlambda_2), dydlambda = 0.125 * F * (sin_2phi * sinlambda_2 - E * sinphi * cos2phi * sinlambda), dydphi = 0.5 * F * (sin2phi * coslambda_2 + E * sin2lambda_2 * cosphi) + 0.5, denominator = dxdphi * dydlambda - dydphi * dxdlambda, dlambda = (fy * dxdphi - fx * dydphi) / denominator, dphi = (fx * dydlambda - fy * dxdlambda) / denominator;
    lambda -= dlambda, phi -= dphi;
  } while ((abs2(dlambda) > epsilon3 || abs2(dphi) > epsilon3) && --i > 0);
  return [lambda, phi];
};
function winkel3_default() {
  return projection(winkel3Raw).scale(158.837);
}

// projections/helpers.ts
function multiplex2(streams) {
  return {
    point(x, y) {
      for (const s of streams)
        s.point(x, y);
    },
    sphere() {
      for (const s of streams)
        s.sphere?.() ?? null;
    },
    lineStart() {
      for (const s of streams)
        s.lineStart();
    },
    lineEnd() {
      for (const s of streams)
        s.lineEnd();
    },
    polygonStart() {
      for (const s of streams)
        s.polygonStart();
    },
    polygonEnd() {
      for (const s of streams)
        s.polygonEnd();
    }
  };
}
function isNil(n) {
  return n === null || n === void 0;
}

// projections/albers-usa-pr.ts
function geoAlbersUsaPr() {
  const epsilon4 = 1e-6;
  let cache, cacheStream, lower48Point, alaskaPoint, hawaiiPoint, puertoRicoPoint, point;
  const lower48 = albers_default(), alaska = conicEqualArea_default().rotate([154, 0]).center([-2, 58.5]).parallels([55, 65]), hawaii = conicEqualArea_default().rotate([157, 0]).center([-3, 19.9]).parallels([8, 18]), puertoRico = conicEqualArea_default().rotate([66, 0]).center([0, 18]).parallels([8, 18]), pointStream = { point: function(x, y) {
    point = [x, y];
  } };
  function albersUsaPr(coordinates) {
    const x = coordinates[0];
    const y = coordinates[1];
    return point = null, (lower48Point.point(x, y), point) || (alaskaPoint.point(x, y), point) || (hawaiiPoint.point(x, y), point) || (puertoRicoPoint.point(x, y), point);
  }
  function reset() {
    cache = cacheStream = null;
    return albersUsaPr;
  }
  albersUsaPr.invert = function(coordinates) {
    const k2 = lower48.scale(), t = lower48.translate(), x = (coordinates[0] - t[0]) / k2, y = (coordinates[1] - t[1]) / k2;
    return (y >= 0.12 && y < 0.234 && x >= -0.425 && x < -0.214 ? alaska : y >= 0.166 && y < 0.234 && x >= -0.214 && x < -0.115 ? hawaii : y >= 0.204 && y < 0.234 && x >= 0.32 && x < 0.38 ? puertoRico : lower48).invert?.(coordinates) ?? null;
  };
  albersUsaPr.stream = function(stream) {
    return cache && cacheStream === stream ? cache : cache = multiplex2([
      lower48.stream(cacheStream = stream),
      alaska.stream(stream),
      hawaii.stream(stream),
      puertoRico.stream(stream)
    ]);
  };
  function precision(precision2) {
    if (isNil(precision2))
      return lower48.precision();
    lower48.precision(precision2), alaska.precision(precision2), hawaii.precision(precision2), puertoRico.precision(precision2);
    return reset();
  }
  function scale(scale2) {
    if (isNil(scale2))
      return lower48.scale();
    lower48.scale(scale2), alaska.scale(scale2 * 0.35), hawaii.scale(scale2), puertoRico.scale(scale2);
    return albersUsaPr.translate(lower48.translate());
  }
  function translate(translate2) {
    if (isNil(translate2))
      return lower48.translate();
    const k2 = lower48.scale(), x = +translate2[0], y = +translate2[1];
    lower48Point = lower48.translate(translate2).clipExtent([[x - 0.455 * k2, y - 0.238 * k2], [x + 0.455 * k2, y + 0.238 * k2]]).stream(pointStream);
    alaskaPoint = alaska.translate([x - 0.307 * k2, y + 0.201 * k2]).clipExtent([[x - 0.425 * k2 + epsilon4, y + 0.12 * k2 + epsilon4], [x - 0.214 * k2 - epsilon4, y + 0.234 * k2 - epsilon4]]).stream(pointStream);
    hawaiiPoint = hawaii.translate([x - 0.205 * k2, y + 0.212 * k2]).clipExtent([[x - 0.214 * k2 + epsilon4, y + 0.166 * k2 + epsilon4], [x - 0.115 * k2 - epsilon4, y + 0.234 * k2 - epsilon4]]).stream(pointStream);
    puertoRicoPoint = puertoRico.translate([x + 0.35 * k2, y + 0.224 * k2]).clipExtent([[x + 0.32 * k2, y + 0.204 * k2], [x + 0.38 * k2, y + 0.234 * k2]]).stream(pointStream);
    return reset();
  }
  albersUsaPr.precision = precision;
  albersUsaPr.scale = scale;
  albersUsaPr.translate = translate;
  return albersUsaPr.scale(1070);
}
var albers_usa_pr_default = geoAlbersUsaPr;

// projections/albers-usa-territories.ts
function geoAlbersUsaTerritories() {
  const epsilon4 = 1e-6;
  let cache, cacheStream, lower48Point, alaskaPoint, hawaiiPoint, puertoRicoPoint, guamMarianaPoint, americanSamoaPoint, point;
  const lower48 = albers_default(), alaska = conicEqualArea_default().rotate([154, 0]).center([-2, 58.5]).parallels([55, 65]), hawaii = conicEqualArea_default().rotate([157, 0]).center([-3, 19.9]).parallels([8, 18]), puertoRico = conicEqualArea_default().rotate([66, 0]).center([0, 18]).parallels([8, 18]), guamMariana = conicEqualArea_default().rotate([-145, 0]).center([0, 16]).parallels([10, 20]), americanSamoa = conicEqualArea_default().rotate([170, 0]).center([0, -14]).parallels([-14, 0]), pointStream = { point: function(x, y) {
    point = [x, y];
  } };
  function albersUsaTerritories(coordinates) {
    const x = coordinates[0];
    const y = coordinates[1];
    point = null;
    return (lower48Point.point(x, y), point) || (alaskaPoint.point(x, y), point) || (hawaiiPoint.point(x, y), point) || (puertoRicoPoint.point(x, y), point) || (guamMarianaPoint.point(x, y), point) || (americanSamoaPoint.point(x, y), point);
  }
  function reset() {
    cache = cacheStream = null;
    return albersUsaTerritories;
  }
  albersUsaTerritories.invert = function(coordinates) {
    const k2 = lower48.scale(), t = lower48.translate(), x = (coordinates[0] - t[0]) / k2, y = (coordinates[1] - t[1]) / k2;
    return (y >= 0.12 && y < 0.234 && x >= -0.39 && x < -0.185 ? alaska : y >= 0.166 && y < 0.234 && x >= -0.185 && x < -0.08 ? hawaii : y >= 0.204 && y < 0.234 && x >= 0.3 && x < 0.38 ? puertoRico : y >= 0.05 && y < 0.21 && x >= -0.45 && x < -0.39 ? guamMariana : y >= 0.21 && y < 0.234 && x >= -0.45 && x < -0.39 ? americanSamoa : lower48).invert?.(coordinates) ?? null;
  };
  albersUsaTerritories.stream = function(stream) {
    return cache && cacheStream === stream ? cache : cache = multiplex2([
      lower48.stream(cacheStream = stream),
      alaska.stream(stream),
      hawaii.stream(stream),
      puertoRico.stream(stream),
      guamMariana.stream(stream),
      americanSamoa.stream(stream)
    ]);
  };
  function precision(precision2) {
    if (isNil(precision2))
      return lower48.precision();
    lower48.precision(precision2);
    alaska.precision(precision2);
    hawaii.precision(precision2);
    puertoRico.precision(precision2);
    guamMariana.precision(precision2);
    americanSamoa.precision(precision2);
    return reset();
  }
  function scale(scale2) {
    if (isNil(scale2))
      return lower48.scale();
    lower48.scale(scale2);
    alaska.scale(scale2 * 0.35);
    hawaii.scale(scale2);
    puertoRico.scale(scale2);
    guamMariana.scale(scale2);
    americanSamoa.scale(scale2);
    return albersUsaTerritories.translate(lower48.translate());
  }
  function translate(translate2) {
    if (isNil(translate2))
      return lower48.translate();
    const k2 = lower48.scale();
    const x = +translate2[0], y = +translate2[1];
    lower48Point = lower48.translate(translate2).clipExtent([[x - 0.455 * k2, y - 0.238 * k2], [x + 0.455 * k2, y + 0.238 * k2]]).stream(pointStream);
    alaskaPoint = alaska.translate([x - 0.275 * k2, y + 0.201 * k2]).clipExtent([[x - 0.39 * k2 + epsilon4, y + 0.12 * k2 + epsilon4], [x - 0.185 * k2 - epsilon4, y + 0.234 * k2 - epsilon4]]).stream(pointStream);
    hawaiiPoint = hawaii.translate([x - 0.18 * k2, y + 0.212 * k2]).clipExtent([[x - 0.185 * k2 + epsilon4, y + 0.166 * k2 + epsilon4], [x - 0.08 * k2 - epsilon4, y + 0.234 * k2 - epsilon4]]).stream(pointStream);
    puertoRicoPoint = puertoRico.translate([x + 0.335 * k2, y + 0.224 * k2]).clipExtent([[x + 0.3 * k2, y + 0.204 * k2], [x + 0.38 * k2, y + 0.234 * k2]]).stream(pointStream);
    guamMarianaPoint = guamMariana.translate([x - 0.415 * k2, y + 0.14 * k2]).clipExtent([[x - 0.45 * k2, y + 0.05 * k2], [x - 0.39 * k2, y + 0.21 * k2]]).stream(pointStream);
    americanSamoaPoint = americanSamoa.translate([x - 0.415 * k2, y + 0.215 * k2]).clipExtent([[x - 0.45 * k2, y + 0.21 * k2], [x - 0.39 * k2, y + 0.234 * k2]]).stream(pointStream);
    return reset();
  }
  albersUsaTerritories.precision = precision;
  albersUsaTerritories.scale = scale;
  albersUsaTerritories.translate = translate;
  return albersUsaTerritories.scale(1070);
}
var albers_usa_territories_default = geoAlbersUsaTerritories;

// projections/index.ts
var R = 6378137;
var projections = {};
var geoProjs = {
  albers: albers_default,
  albersUsa: albersUsa_default,
  azimuthalEqualArea: azimuthalEqualArea_default,
  azimuthalEquidistant: azimuthalEquidistant_default,
  conicConformal: conicConformal_default,
  conicEqualArea: conicEqualArea_default,
  conicEquidistant: conicEquidistant_default,
  equalEarth: equalEarth_default,
  gnomonic: gnomonic_default,
  mercator: mercator_default,
  naturalEarth1: naturalEarth1_default,
  orthographic: orthographic_default,
  stereographic: stereographic_default,
  transverseMercator: transverseMercator_default
};
var geoProjProjs = {
  airy: airy_default,
  aitoff: aitoff_default,
  armadillo: armadillo_default,
  august: august_default,
  baker: baker_default,
  berghaus: berghaus_default,
  // bertin1953: geoBertin1953,
  boggs: boggs_default,
  bonne: bonne_default,
  bottomley: bottomley_default,
  bromley: bromley_default,
  // chamberlin: geoChamberlin,
  collignon: collignon_default,
  craig: craig_default,
  craster: craster_default,
  cylindricalEqualArea: cylindricalEqualArea_default,
  cylindricalStereographic: cylindricalStereographic_default,
  eckert1: eckert1_default,
  eckert2: eckert2_default,
  eckert3: eckert3_default,
  eckert4: eckert4_default,
  eckert5: eckert5_default,
  eckert6: eckert6_default,
  eisenlohr: eisenlohr_default,
  fahey: fahey_default,
  foucaut: foucaut_default,
  // foucautSinusoidal: geoFoucautSinusoidal,
  gilbert: gilbert_default,
  gingery: gingery_default,
  ginzburg4: ginzburg4_default,
  ginzburg5: ginzburg5_default,
  ginzburg6: ginzburg6_default,
  ginzburg8: ginzburg8_default,
  ginzburg9: ginzburg9_default,
  gringorten: gringorten_default,
  guyou: guyou_default,
  hammer: hammer_default,
  hammerRetroazimuthal: hammerRetroazimuthal_default,
  healpix: healpix_default,
  hill: hill_default,
  homolosine: homolosine_default,
  // hufnagel: geoHufnagel,
  // hyperelliptical: geoHyperelliptical,
  // interrupt: geoInterrupt,
  interruptedBoggs: boggs_default2,
  interruptedHomolosine: homolosine_default2,
  interruptedMollweide: mollweide_default2,
  interruptedMollweideHemispheres: mollweideHemispheres_default,
  // interruptedQuarticAuthalic: geoInterruptedQuarticAuthalic,
  interruptedSinuMollweide: sinuMollweide_default2,
  interruptedSinusoidal: sinusoidal_default2,
  kavrayskiy7: kavrayskiy7_default,
  lagrange: lagrange_default,
  larrivee: larrivee_default,
  laskowski: laskowski_default,
  littrow: littrow_default,
  loximuthal: loximuthal_default,
  miller: miller_default,
  // modifiedStereographic: geoModifiedStereographic,
  mollweide: mollweide_default,
  mtFlatPolarParabolic: mtFlatPolarParabolic_default,
  mtFlatPolarQuartic: mtFlatPolarQuartic_default,
  mtFlatPolarSinusoidal: mtFlatPolarSinusoidal_default,
  // naturalEarth2: geoNaturalEarth2,
  nellHammer: nellHammer_default,
  // nicolosi: geoNicolosi,
  patterson: patterson_default,
  polyconic: polyconic_default,
  // polyhedral: geoPolyhedral,
  polyhedralButterfly: butterfly_default,
  polyhedralCollignon: collignon_default2,
  polyhedralWaterman: waterman_default,
  gringortenQuincuncial: gringorten_default2,
  peirceQuincuncial: peirce_default,
  // quincuncial: geoQuincuncial,
  rectangularPolyconic: rectangularPolyconic_default,
  robinson: robinson_default,
  satellite: satellite_default,
  sinuMollweide: sinuMollweide_default,
  sinusoidal: sinusoidal_default,
  times: times_default,
  // twoPointAzimuthal: geoTwoPointAzimuthal,
  // twoPointEquidistant: geoTwoPointEquidistant,
  vanDerGrinten: vanDerGrinten_default,
  vanDerGrinten2: vanDerGrinten2_default,
  vanDerGrinten3: vanDerGrinten3_default,
  vanDerGrinten4: vanDerGrinten4_default,
  // wagner: geoWagner,
  wagner4: wagner4_default,
  wagner6: wagner6_default,
  wiechel: wiechel_default,
  winkel3: winkel3_default
};
var ourProjs = {
  albersUsaPr: albers_usa_pr_default,
  albersUsaTerritories: albers_usa_territories_default
};
Object.keys(geoProjs).forEach((proj) => {
  projections[proj] = geoProjs[proj]().translate([0, 0]).scale(R);
});
Object.keys(geoProjProjs).forEach((proj) => {
  projections[proj] = geoProjProjs[proj]().translate([0, 0]).scale(R);
});
Object.keys(ourProjs).forEach((proj) => {
  projections[proj] = ourProjs[proj]().translate([0, 0]).scale(R);
});
var projections_default = projections;

// src/index.ts
function reproject(options, geometry) {
  const streams = [];
  if (options.forward) {
    let proj = options.forward;
    if (typeof proj === "string") {
      proj = options.projections[proj];
    }
    streams.push(proj.stream, flipY());
  }
  if (options.reverse) {
    let proj = options.reverse;
    if (typeof proj === "string") {
      proj = options.projections[proj];
    }
    streams.push(reverse(proj), flipY());
  }
  streams.reverse();
  const projection2 = {
    stream: function(output) {
      return streams.reduce((combined, s) => s(combined), output);
    }
  };
  return project_default(geometry, projection2);
}
function reverse(projection2) {
  let prev = [];
  return transform_default({
    point: function(x, y) {
      x = clamp(x, -2003750834278924e-8, 2003750834278924e-8);
      const reversed = projection2.invert?.([x, y]);
      if (!reversed || reversed[0] === prev[0] && reversed[1] === prev[1]) {
        return;
      }
      prev = reversed;
      this.stream.point(reversed[0], reversed[1]);
    }
  }).stream;
}
function flipY() {
  return transform_default({
    point: function(x, y) {
      this.stream.point(x, -y);
    }
  }).stream;
}
function clamp(value, min2, max2) {
  if (value > max2) {
    return max2;
  } else if (value < min2) {
    return min2;
  } else {
    return value;
  }
}
export {
  projections_default as projections,
  reproject
};
