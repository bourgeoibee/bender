use std::{ops::{Add, Mul, Sub, Div, Neg}, fmt::{Display, Formatter} };

const T_MAX: f64 = f64::MAX;
const T_MIN: f64 = 0.;

fn main() {
    const ORIGIN: Point = Vec3::new(0., 0., 0.);

    const ASPECT_RATIO: f64 = 16. / 9.;

    const IMAGE_HEIGHT: u64 = 512;
    const IMAGE_WIDTH: u64 = (ASPECT_RATIO * IMAGE_HEIGHT as f64) as u64;

    const FOCAL_LENGTH: f64 = 1.;
    const VIEWPORT_HEIGHT: f64 = 2.;
    const VIEWPORT_WIDTH: f64 = ASPECT_RATIO * VIEWPORT_HEIGHT;

    const VIEWPORT_BOTTOM_LEFT_CORNER: Point = Point::new(-0.5 * VIEWPORT_WIDTH, -0.5 * VIEWPORT_HEIGHT, -FOCAL_LENGTH);
    const HORIZONTAL: Vec3 = Vec3::new(VIEWPORT_WIDTH, 0., 0.);
    const VERTICAL: Vec3 = Vec3::new(0., VIEWPORT_HEIGHT, 0.);

    let scene = vec![
        Shape::Sphere { center: Point::new(0., 0., -2. * FOCAL_LENGTH), radius: 5.},
        //Shape::Sphere { center: Point::new(0., -1., -5. * FOCAL_LENGTH), radius: 5.}
    ];

    println!("P3");
    println!("{IMAGE_WIDTH} {IMAGE_HEIGHT}");
    println!("255");

    for y in (0..IMAGE_HEIGHT).rev() {
        eprintln!("Lines remaining: {}", y);

        let v = (y as f64) / IMAGE_HEIGHT as f64;

        for x in 0..IMAGE_WIDTH {
            let u = (x as f64) / IMAGE_WIDTH as f64;
            let ray = Ray { origin: ORIGIN, direction: VIEWPORT_BOTTOM_LEFT_CORNER + HORIZONTAL * u + VERTICAL * v };
            //dbg!(&ray);
            //std::process::exit(0);
            let color = ray_color(ray, &scene);

            print_color(color);
        }
    }
}

fn ray_color(ray: Ray, scene: &[Shape]) -> Color {
    let mut nearest_t = T_MAX;
    let mut hit_normal: Option<Vec3> = None;

    for shape in scene {
        if let Some(hit_record) = shape.hit(ray, T_MIN, nearest_t) {
            assert!(hit_record.t < T_MIN || hit_record.t > nearest_t);

            nearest_t = hit_record.t;
            hit_normal = Some(hit_record.normal);
        }
    }

    // TODO: make into sky
    match hit_normal {
        Some(normal) => {
            normal * 0.5 + 0.5
        },
        None => {
            const SKY_LOW: Color = Color::new(0.5, 0.9, 0.95);
            const SKY_HIGH: Color = Color::new(0.1, 0.2, 0.6);

            let t = ray.direction.normalized().y * 0.5 + 0.5;

            let sky_color = SKY_LOW * (1. - t) + SKY_HIGH * t;
            //dbg!(sky_color);
            sky_color
        }
    }
}

fn dot(lhs: Vec3, rhs: Vec3) -> f64 {
    lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z
}

fn print_color(color: Color) {
    let color = color * 255.;
    let red = color.x as u8;
    let green = color.y as u8;
    let blue = color.z as u8;

    println!("{red} {green} {blue}");
}

#[derive(Debug, Clone, Copy)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

struct HitRecord {
    point: Point,
    normal: Vec3,
    t: f64,
}

struct Scene(Vec<Box<dyn SceneObject>>);

type Point = Vec3;
type Color = Vec3;

trait SceneObject {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

impl Vec3 {
    const fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    fn squared_length(self) -> f64 {
        dot(self, self)
    }

    fn normalized(self) -> Self {
        self / self.squared_length().sqrt()
    }
}

impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, rhs: Self) -> Self {
        Vec3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Add<f64> for Vec3 {
    type Output = Vec3;

    fn add(self, rhs: f64) -> Self {
        Vec3::new(self.x + rhs, self.y + rhs, self.z + rhs)
    }
}

impl Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, rhs: Self) -> Self {
        Vec3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Sub<f64> for Vec3 {
    type Output = Vec3;

    fn sub(self, rhs: f64) -> Self {
        Vec3::new(self.x - rhs, self.y - rhs, self.z - rhs)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;

    fn mul(self, rhs: f64) -> Self {
        Vec3::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;

    fn div(self, rhs: f64) -> Self {
        Vec3::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Self {
        self * -1.
    }
}

#[derive(Debug, Clone, Copy)]
struct Ray {
    origin: Vec3,
    direction: Vec3,
}

impl Ray {
    fn at(self, t: f64) -> Point {
        self.origin + self.direction * t
    }
}

enum Shape {
    Sphere {
        center: Point,
        radius: f64,
    },
}

impl Shape {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        match self {
            Shape::Sphere { center, radius } => {
                // Uses the quadratic formula
                let circle_to_camera = ray.origin - *center;
                let a = ray.direction.squared_length();
                let half_b = dot(circle_to_camera, ray.direction);
                let c = circle_to_camera.squared_length() - radius * radius;

                let discrminant = half_b * half_b - a * c;
                if discrminant < 0. { return None }
                
                let t = (-half_b - discrminant.sqrt()) / a;
                if t < t_min || t > t_max {
                    let t = (-half_b + discrminant.sqrt()) / a;
                    if t < t_min || t > t_max {
                        return None;
                    }
                 }

                let point = ray.at(t);
                dbg!(point);
                Some(HitRecord { t, point, normal: (point - *center) / *radius })
            },
        }
    }
}

//impl Scene {
//    fn new() -> Self {
//        Scene(vec![])
//    }
//
//    fn push(&mut self, object: &dyn SceneObject) {
//        self.0.push(Box::new(*object));
