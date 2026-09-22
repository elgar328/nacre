//! Fixture profiles stated as edge lists (`line`, `arc_turns`, `arc_rat`, `circle`), **ordered**.
//!
//! Every fixture written with it is a chain: each edge starting where the previous one ended, a
//! circle on its own. The kernel takes rings, not soups (`Ring2d::new`, `from_paths`), so this
//! walks the edges in the order they were written and hands the rings over. Nothing here
//! re-discovers an order; an edge that does not start where the pen stands is a fixture bug and
//! panics. Ring order is whole circles first, then the chains as written, each seeded at its
//! first edge's start.
#![allow(dead_code)]

use nacre_exact::Rat;
use nacre_math::Point2;
use nacre_ops::{Edge2d, Profile2d, Ring2d, SketchError, from_paths};

/// One edge as a fixture spells it.
#[derive(Clone, Copy, Debug)]
pub enum Stated {
    Line(Point2, Point2),
    Arc {
        center: Point2,
        start: Point2,
        turns: i32,
    },
    ArcRat {
        center: [Rat; 2],
        start: [Rat; 2],
        end: [Rat; 2],
        ccw: bool,
    },
    Circle {
        center: Point2,
        radius: f64,
    },
}

pub fn line(from: Point2, to: Point2) -> Stated {
    Stated::Line(from, to)
}

pub fn arc_turns(center: Point2, start: Point2, turns: i32) -> Stated {
    Stated::Arc {
        center,
        start,
        turns,
    }
}

pub fn arc_rat(center: [Rat; 2], start: [Rat; 2], end: [Rat; 2], ccw: bool) -> Stated {
    Stated::ArcRat {
        center,
        start,
        end,
        ccw,
    }
}

pub fn circle(center: Point2, radius: f64) -> Stated {
    Stated::Circle { center, radius }
}

fn lift(p: Point2) -> Result<[Rat; 2], SketchError> {
    match (Rat::from_decimal(p[0]), Rat::from_decimal(p[1])) {
        (Some(x), Some(y)) => Ok([x, y]),
        _ => Err(SketchError::OutsideDecimalWindow { at: p.as_array() }),
    }
}

/// The edges in the order written, walked into rings, then sorted into profiles.
pub fn stated(edges: Vec<Stated>) -> Result<Vec<Profile2d>, SketchError> {
    let mut circles: Vec<Ring2d> = Vec::new();
    let mut chains: Vec<Ring2d> = Vec::new();
    let (mut vertices, mut steps): (Vec<[Rat; 2]>, Vec<Edge2d>) = (Vec::new(), Vec::new());
    let mut here: Option<[Rat; 2]> = None;
    for e in edges {
        let (from, step, to) = match e {
            Stated::Circle { center, radius } => {
                circles.push(Ring2d::circle(center, radius)?);
                continue;
            }
            Stated::Line(a, b) => (lift(a)?, Edge2d::Line, lift(b)?),
            Stated::Arc {
                center,
                start,
                turns,
            } => {
                let (c, s) = (lift(center)?, lift(start)?);
                let (step, end) = nacre_ops::arc_turns_rat(c, s, turns)?;
                (s, step, end)
            }
            Stated::ArcRat {
                center,
                start,
                end,
                ccw,
            } => (start, nacre_ops::arc_to_rat(center, start, end, ccw)?, end),
        };
        if let Some(h) = here {
            assert_eq!(
                h, from,
                "a fixture edge must start where the previous one ended"
            );
        }
        vertices.push(from);
        steps.push(step);
        here = Some(to);
        if to == vertices[0] {
            chains.push(Ring2d::new(
                std::mem::take(&mut vertices),
                std::mem::take(&mut steps),
            )?);
            here = None;
        }
    }
    assert!(here.is_none(), "a fixture chain did not close");
    circles.extend(chains);
    from_paths(circles)
}
