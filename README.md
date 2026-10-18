# interpolate

Linear interpolation for [Tessel](https://tessel-lang.syncster.dev)
programs: given two known points, find the x that goes with a y between
them, the way you read a value off a table.

```tessel
let a = interpolate.Point(x: 10.0, y: 100.0)
let b = interpolate.Point(x: 20.0, y: 200.0)
interpolate.interpolate(p1: a, p3: b, y2: 150.0)   // 15.0
```

## Using it

Add it to a project (Tessel 0.8 or later):

```sh
tessel add interpolate
```

then, in the program:

```tessel
import interpolate

fn main() {
    // A rating table: water depth (x, in m) for two discharges (y, in m³/s).
    let low = interpolate.Point(x: 1.20, y: 35.0)
    let high = interpolate.Point(x: 1.80, y: 62.0)

    // The depth for a discharge of 50 m³/s.
    let depth = interpolate.interpolate(p1: low, p3: high, y2: 50.0)
    print("depth: {depth} m")
}
```

```text
depth: 1.5333333333333332 m
```

## What it has

### `Point`

```tessel
public struct Point {
    x: Float,
    y: Float
}
```

A known point of the relation: an `x` and the `y` that goes with it. Make
one with `interpolate.Point(x: 10.0, y: 100.0)`.

### `interpolate(p1:p3:y2:)`

```tessel
public fn interpolate(p1: Point, p3: Point, y2: Float) -> Float
```

The x of the point on the straight line through `p1` and `p3` whose y is
`y2`:

```text
x2 = (x1 - x3) / (y1 - y3) * (y2 - y3) + x3
```

| Parameter | What it is |
|---|---|
| `p1` | One known point (point 1). |
| `p3` | The other known point (point 3). |
| `y2` | The y of the point to find (point 2), usually between `p1.y` and `p3.y`. |

It returns the x of point 2.

- **Beyond the points:** a `y2` outside the two points extends the line
  (extrapolation): with the points above, `y2: 250.0` gives `25.0`.
- **Either order:** swapping `p1` and `p3` gives the same answer.
- **Level lines:** the points must have different y's. When `p1.y ==
  p3.y`, no single x has height `y2`, and the result is `inf` or `-inf`
  (or `NaN` when `y2` is that same height). Check for it first if your
  data can have equal y's.

## Publishing it

The library is published from its git repository, on GitHub or another
git host. Commit the code, then:

```sh
tessel publish --dry-run   # check it first
tessel publish             # tag the version, push the tag, propose it to the index
```

Before the next version, raise `version` in `tessel.toml`. See
[Libraries](https://tessel-lang.syncster.dev/guides/libraries/) for how
publishing works.

## What's where

| File | |
|---|---|
| `interpolate.tsl` | The library's code: `Point` and `interpolate`. |
| `tessel.toml` | The library's name, version, description and licence. |

## Licence

MIT.
