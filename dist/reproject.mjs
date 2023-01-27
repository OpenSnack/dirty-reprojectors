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
    var pi = ring[i], xi = pi[0], yi = pi[1], pj = ring[j], xj = pj[0], yj = pj[1];
    if (yi > y ^ yj > y && x < (xj - xi) * (y - yi) / (yj - yi) + xi)
      contains = !contains;
  }
  return contains;
}

// node_modules/d3-geo-projection/src/project/index.js
function project_default(object, projection) {
  var stream = projection.stream, project;
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

// index.ts
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
  const projection = {
    stream: function(output) {
      return streams.reduce((combined, s) => s(combined), output);
    }
  };
  return project_default(geometry, projection);
}
function reverse(projection) {
  let prev = [];
  return transform_default({
    point: function(x, y) {
      x = clamp(x, -2003750834278924e-8, 2003750834278924e-8);
      const reversed = projection.invert?.([x, y]);
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
function clamp(value, min, max) {
  if (value > max) {
    return max;
  } else if (value < min) {
    return min;
  } else {
    return value;
  }
}
export {
  reproject
};
