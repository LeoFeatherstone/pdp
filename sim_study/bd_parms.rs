fn bd_parms(p: f32, delta: f32, r_naught: f32) -> Vec<f32> {
    let psi: f32 = (p * delta) / (1.0-p);
    let lambda: f32 = ((delta) / (1.0-p)) * r_naught;
    
    let mut vec: Vec<f32> = Vec::new();
    
    vec.push(psi);
    vec.push(lambda);
    
    return vec
}

fn main() {
    bd_parms(0.8, 365.25 / 10.0, 2.5);
    println!(" SARS-CoV-2 Parms: {:?}", bd_parms(0.8, 365.25 / 10.0, 2.5));
    println!(" TB Parms: {:?}", bd_parms(0.9, 0.5, 3.0));
    println!(" Shigella spp.: {:?}", bd_parms(0.4, 365.25 / 7.0, 2.0));
}