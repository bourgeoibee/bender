use std::ops::{Add, Mul};

fn main() {
    const IMAGE_HEIGHT: u32 = 512;
    const IMAGE_WIDTH: u32 = 512;

    println!("P3");
    println!("{IMAGE_WIDTH} {IMAGE_HEIGHT}");
    println!("255");

    for y in 0..IMAGE_HEIGHT {
        eprintln!("Lines remaining: {}", IMAGE_HEIGHT - y);
        for x in 0..IMAGE_WIDTH {
            let red = 256 * x / IMAGE_WIDTH;
            let green = 256 * y / IMAGE_HEIGHT;
            let blue = 64;

            println!("{red} {green} {blue}");
        }
    }
}

struct Vec3 {
    x: f32,
    y: f32,
    z: f32,
}

struct Point(Vec3);

struct Color(Vec3);

struct Ray {
    origin: Vec3,
    direction: Vec3,
}

impl Vec3 {
    fn new(x: f32, y: f32, z: f32) -> Self {
        Vec3 { x, y, z }
    }
}

impl Add for Vec3 {
    type Output = Vec3;

    fn add(self: Self, rhs: Self) -> Self {
        Vec3 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Mul<f32> for Vec3 {
    type Output = Vec3;

    fn mul(self, rhs: f32) -> Self {
        Vec3 {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Ray {
    fn at(self, time: f32) -> Point {
        Point(self.origin + self.direction * time)
    }
}
