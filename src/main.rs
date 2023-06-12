use std::ops::{Add, Mul, Sub, Div, Neg};

const T_MAX: f64 = f64::MAX;
const T_MIN: f64 = 0.;

fn main() {
    const ORIGIN: Point = Vec3::new(0., 0., 0.);

    const ASPECT_RATIO: f64 = 16. / 9.;

    const IMAGE_HEIGHT: u64 = 512;
    const IMAGE_WIDTH: u64 = (ASPECT_RATIO * IMAGE_HEIGHT as f64) as u64;

    const VIEWPORT_HEIGHT: f64 = 2.;
    const VIEWPORT_WIDTH: f64 = ASPECT_RATIO * VIEWPORT_HEIGHT;
    const FOCAL_LENGTH: f64 = 1.;

    const VIEWPORT_BOTTOM_LEFT_CORNER: Point = Point::new(-0.5 * VIEWPORT_WIDTH, -0.5 * VIEWPORT_HEIGHT, -FOCAL_LENGTH);

    let scene = [ Box::new(Sphere::new(Point::new(0., 0., -FOCAL_LENGTH), 1.)) ];

    println!("P3");
    println!("{IMAGE_WIDTH} {IMAGE_HEIGHT}");
    println!("255");

    for y in 0..IMAGE_HEIGHT {
        eprintln!("Lines remaining: {}", IMAGE_HEIGHT - y);
        for x in 0..IMAGE_WIDTH {
            let ray = Ray { origin: ORIGIN, direction: VIEWPORT_BOTTOM_LEFT_CORNER + u * HORIZONTAL + v * VERTICAL };
            let color = ray_color(ray, scene);

            println!("{}", color);
        }
    }
}

fn ray_color(ray: Ray, scene: &[Box<dyn SceneObject>]) -> Color {
    //let nearest_t = T_MAX;
    let hit: Option<HitRecord> = None;

    for object in scene {
        if let Some(hit_record) = (*object).hit(ray, T_MIN, nearest_t) {
            if hit_record.t < T_MIN || hit_record.t > nearest_t {
                continue;
            }

            hit = Some(hit_record);
        }
    }

    // TODO: make into sky
    match hit {
        Some(record) => {
            let color = (record.normal + 1) * 0.5;
        }
        None => Color::new(0., 0., 0.);
    }


}

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

type Point = Vec3;
type Color = Vec3;

trait SceneObject {
    fn hit(self, ray: Ray, t_min: f64, t_max: f64) -> Option(HitRecord);
}

impl Vec3 {
    const fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    fn squared_length(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
}

impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, rhs: Self) -> Self {
        Vec3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, rhs: Self) -> Self {
        Vec3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
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

impl Display for Vec3 {
    
}

struct Ray {
    origin: Vec3,
    direction: Vec3,
}

impl Ray {
    fn at(self, time: f64) -> Point {
        self.origin + self.direction * time
    }
}

struct Sphere {
    center: Point,
    radius: f64,
}

impl Sphere {
    fn new(center: Point, radius: f64) -> Self {
        Sphere { center, radius }
    }
}

impl SceneObject for Sphere {
    // Uses the quadratic formula
    fn hit(self, ray: Ray, t_min: f32, _max: f32) -> Option(HitRecord) {
        let center_to_tail = ray.origin - self.center;
        let a = ray.direction.squared_length();
        let half_b = dot(center_to_tail, ray.direction);
        let c = center_to_tail.squared_length() - radius * radius;

        let discrminant = half_b * half_b - a * c;
        if discrminant < 0 { return None }
        
        let t = (-half_b - discrminant.sqrt()) / a;
        if t < t_min || t > t_max {
            let t = (-half_b + discrminant.sqrt()) / a;
            if t < t_min || t > t_max {
                return None;
            }
         }

        let point = ray.at(t);
        Some(HitRecord { t, point, normal: (point - self.center) / self.radius })
    }
}

