
/******************************************************************************* 

2021 Dr. Berke et Dr. Ben 

: cluster (tSNE) by curvature (angeldefect) along the main diagonal of
3D Models of pottery classes 

I AM NOT WILLING TO THINK THAT WAY. I STICK WITH NOT ABSTRACT ABSTRACTIONS.

*******************************************************************************/

//use glium::glutin::dpi::LogicalSize;
//use glium::glutin::event::{Event, WindowEvent};
//use glium::glutin::event_loop::{ControlFlow, EventLoop};
//use glium::glutin::window::WindowBuilder;
//use glium::glutin::ContextBuilder;
//use glium::uniform;
//use glium::Program;
use tobj;
use num_cpus;
use std::fs;
use std::env;
use core::f32::consts::PI;
use std::collections::HashMap;
use std::sync::mpsc;
use std::thread;
use rand::Rng;
extern crate serde;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/*---------------------geom datatypes-----------------------------------------*/
#[derive(Debug, PartialEq, PartialOrd)]
struct Tristar {
    poi: [f32; 3],
    sta: Vec<[[f32; 3]; 3]>
}

impl Tristar {
    pub fn new(poi: [f32; 3], sta: Vec<[[f32; 3]; 3]>) -> Self {
        Tristar {
            poi,
            sta
        }
    }
}

/*---------------------vec                   arithmetics----------------------*/
fn p32vec3( p1: &[f32; 3], p2: &[f32; 3] ) -> [f32; 3]{
    return [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2] ];
}

fn vec3len( thevec: &[f32; 3] ) -> f32 {
    return ( thevec[0].powi( 2 ) + thevec[1].powi( 2 ) + thevec[2].powi( 2 ) ).sqrt( );
}

fn vec3minus( v1: &[f32; 3], v2: &[f32; 3] ) -> [f32; 3] {
    return p32vec3(v1, v2);
}

fn vec3divscaler( v: &[f32; 3], s: &f32 ) -> [f32; 3] {
    return [v[0]/s, v[1]/s, v[2]/s];
}

fn vec3crossvec3( v1: &[f32; 3], v2: &[f32; 3] ) -> [f32; 3] {
    return  [v1[1]*v2[0] - v1[2]*v2[1] , v1[0]*v2[0]-v1[2]*v2[2], v1[0]*v2[1]-v1[1]*v2[2]]
}

fn vec3innerprodvec3( v1: &[f32; 3], v2: &[f32; 3] ) -> f32 {
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

fn vec3angelvec3( v1: &[f32; 3], v2: &[f32; 3] ) -> f32 {
    //return ( vec3innerprodvec3( v1, v2 ) / ( vec3len( v1 ) * vec3len( v2 ) ) ).acos( );
    return vec3innerprodvec3( v1, v2 ) / ( vec3len( v1 ) * vec3len( v2 ) );
}

fn vec3eqvec3( v1: &[f32; 3], v2: &[f32; 3] ) -> bool {
    if v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2] {
        return true;
    } else {
        return false;
    }
}

fn compdistp3vec3( p1: &[f32; 3], p2: &[f32; 3], v1: &[f32; 3], l1: &f32 ) -> f32 {
    let to = p32vec3( p2, p1 );
    let cp = vec3crossvec3( &to, &v1 ); 
    let d = vec3len( &cp ) / l1; 
    return d;
}

fn compdistvecnvecn( v1: &[f32], v2: &[f32]) -> f32 {
    //println!("{:?}---{:?}", v1.len(), v2.len() );
    let mut d: f32 = 0.0;
    for i in 0..v1.len( ) {
        let td = v1[i] - v2[i];
        let qtd = td*td;
        d += qtd;
    }
    return d.sqrt();
}

/*----------------------k-form arithmetic-------------------------------------*/
fn angeldefect( astar: &Vec<[[f32; 3]; 3]>, centerp: &[f32; 3] ) -> f32 {
    let mut starang = 0.0;
    let mut allnormals = Vec::new();
    for trii in 0..astar.len() {
        let nelem = [0.0, 0.0, 0.0];
        let mut a = [0.0, 0.0, 0.0];
        let mut b = [0.0, 0.0, 0.0];
        for h in 0..astar[ trii ].len() {
            if !vec3eqvec3( &astar[trii][h], centerp) {
                if vec3eqvec3( &a, &nelem ) {
                    a = astar[ trii ][ h ];
                } else if vec3eqvec3( &b, &nelem ) {
                    b = astar[ trii ][ h ];
                } else {
                    println!("nooooooooo sth wrong")
                }
            }
        }
        //compute normal on dreieck, take points as ortsvectors
        let nz1 = vec3minus( &a, &centerp );
        let nz2 = vec3minus( &b, &centerp );
        let no = vec3crossvec3( &nz1, &nz2 );
        allnormals.push( no );
        let e1 = p32vec3( &centerp, &a );
        let e2 = p32vec3( &centerp, &b );
        starang += vec3angelvec3( &e1, &e2 );
    }

    //ceck if normals will intersect - not concex
    /*let mut alldifflen = 0.0;
    let mut alllen = 0.0;
    for f in 1..allnormals.len() {
        alldifflen += vec3len(&vec3minus(&allnormals[f-1], &allnormals[f]));
        alllen += vec3len(&allnormals[f-1]);
    }
    alllen += vec3len(&allnormals[allnormals.len()-1]);
    let mut orientationadd = 1.0;
    if alldifflen > alllen {
        orientationadd = -1.0;
    }
    let angdef = (2.0*PI - starang)*orientationadd;*/
    let angdef = 2.0*PI - starang;
    return angdef; //is not signed no notion of direction, good or bad????
}

/*---------------------ARRAY / LIST / VEC Type minmax                 --------*/
fn getarrayofnimwithval( n:usize, m: usize, v: f32 ) -> Vec<Vec<f32>> {
    let mut retarr = Vec::new();
    for _i in 0..n {
        let mut retvec = Vec::new();
        for _j in 0 ..m {
            retvec.push( v ); //clone???
        }
        retarr.push( retvec );
    }
    return retarr; 
}

fn zeros( dim: usize ) -> Vec<f32> {
    let mut rv = Vec::new();
    for _i in 0..dim {
        rv.push( 0.0 );
    }
    return rv;
}

fn zeroszeros( n: usize, dim: usize ) -> Vec<Vec<f32>> {
    return getarrayofnimwithval( n, dim, 0.0 );
}

fn iteminarray( arr: &Vec<[f32; 3]>, val: &[f32; 3] ) -> bool {

    for i in 0.. arr.len() {
        if val[0] == arr[i][0] &&  val[1] == arr[i][1]  && val[2] == arr[i][2]  {
            return true;
        }
    }
    return false;
}

fn getmaxdiffpos( a: &[f32] ) -> (usize, f32) {
    let mut absmaxindex = 0; 
    let mut absmaxat = 0.0;
    let mut taken = 0.0;
    for i in 1..a.len() {
        let t = a[i] - a[i-1];
        if absmaxat < t {
            absmaxat = t;
            taken = a[i-1];
            absmaxindex = i;
        }
        
    } 
    return (absmaxindex, taken);
}

fn getminvalandindex( a: &[f32] ) -> (usize, f32) {
    let mut minindex = 0; 
    let mut absmin:f32 = 0.0;
    for i in 0..a.len() {
        if absmin.abs() > a[i].abs() { //minimal numbernot minimal of all 
            absmin = a[i];
            minindex = i;
        }
        
    } 
    return (minindex, absmin);
}


/* advanced distance computation **********************************************/
fn dmultidimvecn( v1: &Vec<f32>, v2: &Vec<f32> ) -> f32 {
    //vergleiche zwei vectoren und mach diese an ihrer mitte fest
    let mut vv1 = v1;
    let mut vv2 = v2;
    
    if v1.len() > v2.len() { 
        vv1 = v2;
        vv2 = v1;
    }

    let vv1l = vv1.len( );
    let offset = (( vv2.len( ) - vv1l ) / 2) as usize;
    return compdistvecnvecn( &vv1[..], &vv2[ offset..(vv1l+offset) ] );
}

fn wasserst1d( v1: &Vec<f32>, v2: &Vec<f32> ) -> f32 {
    //vergleiche zwei vectoren und mach diese an ihrer mitte fest
    let mut vv1 = v1;
    let mut vv2 = v2;
    
    if v1.len() < v2.len() { 
        vv1 = v2;
        vv2 = v1;
    }

    let mut lendiff = vv1.len() - vv2.len();
    let mut intsertafter = 0.0;
    let n = vv1.len( );
    let mut m = vv2.len( );
    if n == 0 {
        return std::f32::INFINITY;
    }
    if m == 0 {
        return std::f32::INFINITY;
    }
    if lendiff != 0 {
        intsertafter = m as f32 / lendiff as f32; //usize got an implicit floor() ????
    }
    //println!("intsertafter: {:?}, n: {:?}, m: {:?}, lendiff, {:?}", intsertafter, n, m, lendiff);
    if lendiff != 0 && intsertafter == 0.0 {
        return std::f32::INFINITY;
    }

    //DO INSERTIONS EQUALY DISTRIBUTED ON THE SMALLER ARRAY
    let mut vv2filled = Vec::new( );
    if intsertafter < 1.0 &&  intsertafter > 0.0 {
        let howmuchfill = (1.0/intsertafter).round() as usize;
        for j in 0..m {
            vv2filled.push( vv2[j] );
            for p in 0..howmuchfill {
                vv2filled.push( vv2[j] );
            }
        }
        if vv2filled.len() < n {
            lendiff = n - vv2filled.len();
            for o in 0..lendiff {
                vv2filled.push(vv2[m-1]);
            }
        }
    }
    if intsertafter >= 1.0 {
        let insertstep = intsertafter.floor() as usize;
        for j in 0..m {
            vv2filled.push( vv2[j] );
            if j % insertstep == 0 {
                vv2filled.push( vv2[j] );
            }
        }
    }
    m = vv2filled.len();
    let mut transportcost:f32 = 0.0;
    let mut nichtverbracht:f32 = 0.0;
    //println!(" m:, {:?} n: {:?}", m, n);
    let mut sizetoitterover = m; //this is needed because of the rounding of the insertion step
    if( m > n ){
        sizetoitterover = n;
    }
    for j in 0..sizetoitterover {
        let mut zeile = Vec::new( );
        for i in 0..sizetoitterover { //should be n, but n should be m sized now
            let u = vv2filled[j] - 
            vv1[i];//unterschied der krümmungen, mit Vorzeichen!!!
            zeile.push( u + nichtverbracht );
           
        }
        //ort des minimalen unterschiedes zur aktuellen position im Array (j zeilennummer)
        let (minvalindex, minval) = getminvalandindex( &zeile );
        let ortdiff = minval - (j as f32);
        //minimaler Unterschied zu nächstem unterschied 
        
        let valatthis = zeile[j];
        //quasi geornete (keine Sprünge) Transportkosten aufadieren
        transportcost += sigmo((minval-valatthis).abs()) * ortdiff;
        //nicht verbrachte differenz
        nichtverbracht = valatthis;
    }
    //nicht verachten rest mit einrechnen
    transportcost = transportcost + nichtverbracht;
    return transportcost;
}


fn distmultidimvecarray( comp: &Vec<Vec<f32>> ) -> Vec<Vec<f32>> {
    let mut dij = Vec::new( );//Vec<Vec<f32>>; size of comp
    let n = comp.len();
    for i in 0..n {
        let mut di = Vec::new( );
        for j in 0..n {
            if i != j {
                println!("Längen der Vektoren {:?}, {:?}, {:?}, {:?}", i, comp[i].len(), j , comp[j].len() );
                di.push( wasserst1d( &comp[i], &comp[j] ) );
                //di.push( dmultidimvecn( &comp[i], &comp[j] ) );
            } else {
                di.push( 0.0 );
            }
        }
        /*susu = susu/ (n as f32);
        for j in 0..n {
            if i != j {
                di[j] = di[j]/susu;
            }

        }*/
        dij.push( di );
    } 
    return dij;
}

/*---------------------MATH helper                                    --------*/


fn sigmo( val: f32 ) -> f32 {
    const g: f32 = 6.0;

    return 1.0 / ( 1.0 + ( -val / g).exp( ) );

}

fn gausssrand( rng: &mut rand::rngs::ThreadRng ) -> f32 {
    let y: f32 = rng.gen();
    let u: f32 = 2.0 * y - 1.0;
    let x: f32 = rng.gen();
    let v: f32 = 2.0 * x - 1.0;
    let r: f32 = u*u + v*v;
    if r == 0.0 || r > 1.0 { 
        return gausssrand( rng );
    }
    let c = (-2.0*r.log(10.0)/r).sqrt();
    return (u*c+ v*c)/2.0;
}

fn raxxxxxxxndn( mu:f32, std:f32, rng: &mut rand::rngs::ThreadRng ) -> f32 {
    return mu+gausssrand(rng)*std;
}

fn getrandnorm( n:usize, m:usize, rng: &mut rand::rngs::ThreadRng ) -> Vec<Vec<f32>> {
    let mut retarr = Vec::new();
    for _i in 0..n {
        let mut retvec = Vec::new();
        for _j in 0 ..m {
            retvec.push( raxxxxxxxndn(0.0, 0.0001, rng) ); 
        }
        retarr.push( retvec );
    }
    return retarr; 
}

fn sign( v: f32 ) -> usize {
    if v < 0.0  {
        return 0; //-1 return
    } 
    return 1;
}


/*---------------------                      derived space, embeddings--------*/
fn d2p( darr: &Vec<Vec<f32>>, howmanyneig: usize, maxtries: usize, tolller: f32) -> Vec<Vec<f32>> {
    let Htarget = (howmanyneig as f32).log(10.0); // target entropy of distribution
    let n = darr.len();
    let mut pij = Vec::new( ); //replace with array of fixed length

    for i in 0..n {
        let mut perrow = zeros(n); //replace with array of fixed length
        
        let mut notdone = true;
        let mut betamin = -std::f32::INFINITY;
        let mut betamax = std::f32::INFINITY;
        let mut beta = 1.0; // initial value of precision
        let mut num = 1;
        while notdone  {
            //compute bedingte wahrscheinlichkeit
            let mut psum: f32 = 0.0;
            for j in 0..n {
                let mut pt: f32 = 0.0; 
                if darr[i].len() == 0 {
                    pt = 0.0;
                } else if  i != j { 
                    //println!("{:?}",darr[i][j]);
                    pt = (-1.0 * darr[i][j] * beta ).exp();
                }
                perrow[j] = pt;
                psum += pt;
            }
            //normalize pij
            let mut Hhere = 0.0;
            for j in 0..n {
                let mut pt = 0.0;
                if psum != 0.0 {
                    pt = perrow[j] / psum;
                } 
                perrow[j] = pt;
                if pt > 0.0000001 { 
                    Hhere -= pt * pt.log(10.0);
                }
            }
            //check if distribution is peaky or too peaky
            if Hhere > Htarget {
                //entropy was too high (distribution too diffuse)
                //increase the beta (precision) for more peaky distribution
                betamin = beta; // move up the bounds
                if betamax == std::f32::INFINITY { 
                    beta = beta * 2.0; 
                } else { 
                    beta = (beta + betamax) / 2.0; 
                }
            } else {
                //make distrubtion less peaky
                betamax = beta;
                if betamin == -std::f32::INFINITY { 
                    beta = beta / 2.0; 
                } else { 
                    beta = (beta + betamin) / 2.0; 
                }
            }            
            //println!("h: {:?}, beta: {:?} htar {:?}", Hhere, beta, Htarget);
            //quit while
            if (Hhere - Htarget).abs() < tolller { 
                notdone = false; 
            }
            if num >= maxtries { 
                notdone = false; 
            }
            num += 1;
        }
        //println!("num: {:?}", num);   
            
        pij.push( perrow )
    }
    return pij;
}

fn costandgradient( n: usize, dim: usize, step: usize, pij: &Vec<Vec<f32>>, solY: &mut Vec<Vec<f32>> ) -> (f32, Vec<Vec<f32>>) {
    let mut pmul = 1.0;
    if step < 100 {
        pmul = 4.0;
    }
    let mut Qu = zeroszeros( n, n );
    let mut qsum = 0.0;
    for i in 0..n {
        for j in 0..n {
            if Qu[i][j] == 0.0 { // because you fill in both ij and ji in Qu, 
                let mut dsum = 0.0;
                for d in 0..dim {
                    let dhere = solY[i][d] - solY[j][d];
                    dsum += dhere * dhere;
                }
                let qu = 1.0 / (1.0 + dsum); // Student t-distribution
                Qu[i][j] = qu;
                Qu[j][i] = qu;
                qsum += 2.0 * qu;
            }
        }
    }
    //normalize
    let mut Q = zeroszeros( n, n );
    for i in 0..n {
        for j in 0..n {
            Q[i][j] = Qu[i][j] / qsum;
        }
    }
    
    let mut cost = 0.0;
    let mut grad = Vec::new();
    for i in 0..n {
        let mut gsum = zeros( dim );
        for j in 0..n {
            cost += - pij[i][j] * Q[i][j].log(10.0); // accumulate cost
            let premult = 4.0 * (pmul * pij[i][j] - Q[i][j]) * Qu[i][j];
            for d in 0..dim {
                gsum[d] += premult * (solY[i][d] - solY[j][d]);
            }
        }
        grad.push( gsum );
    }
    return (cost, grad);
}



fn optistep( n: usize, dim: usize, step: usize, epsilon: f32, solY: &mut Vec<Vec<f32>>, ystep: &mut Vec<Vec<f32>>, gains: &mut Vec<Vec<f32>>, pij:  &Vec<Vec<f32>>, lc: f32 ) -> f32 {
    //n: usize, dim: usize, step: usize, pij: &Vec<Vec<f32>>, Y: &Vec<Vec<f32>>

    let (cost, grad) = costandgradient( n, dim, step, pij, solY );
    /*if lc < cost {
        return lc;
    }*/
    let mut ymean = zeros( dim );
    for i in 0..n {
        for d in 0..dim {
            let gid = grad[i][d];
            let sid = ystep[i][d];
            let gainid = gains[i][d];

            // compute gain update
            let mut newgain = gainid + 0.2;
            if sign(gid) == sign(sid) {
                 newgain = gainid * 0.8;
            }
            if newgain < 0.01 { 
                newgain = 0.01; // clamp
            }
            gains[i][d] = newgain; // store for next turn

            // compute momentum step diection
            let mut momval = 0.8;
            if step < 250 { 
                momval = 0.5;
            }
            let newsid = momval * sid - epsilon * newgain * grad[i][d];
            ystep[i][d] = newsid; // remember the step we took

            solY[i][d] += newsid; 
            ymean[d] += solY[i][d]; // accumulate mean so that we can center later
        }
    }
    // reproject Y to be zero mean
    for i in 0..n {
        for d in 0..dim {
            solY[i][d] -= ymean[d] / (n as f32);
        }
    }
    return cost;
}

fn tsne( comp: &Vec<Vec<f32>>, dim: usize, pris: usize ) -> Vec<Vec<f32>> {
    //Wann wart ihr jemals dankbar für eure warmen Ärsche im Winter?
    let dij = distmultidimvecarray( comp );
    println!( "dij: {:?}", dij );
    let pij = d2p( &dij, pris, 1000, 0.1 ); //d matrix, neighbors, maxtries of finding good fitt of distribution, tollerance
    println!( "pij: {:?}", pij );
    let n = comp.len();
    let mut rng: rand::rngs::ThreadRng = rand::thread_rng();
    let mut solY = getrandnorm( n, dim, &mut rng );
    println!( "Init Y: {:?}", solY );
    let mut gain = getarrayofnimwithval( n, dim, 1.0 ); // step gains to accelerate progress in unchanging directions
    println!( "Init gain: {:?}", gain );
    let mut ystep = getarrayofnimwithval( n, dim, 0.0 ); //momentum accelerator i.e. step last step memeory
    println!( "Init ystep: {:?}", ystep );

    let mut opticount = 0;
    let maxopti = 1000;
    let mut keeponopti = true;
    let epsilon = 5.0; //"learning" rate 
    let mut lastcost = 1000000.0;
    while keeponopti {
        //n: usize, dim: usize, step: usize, epsilon: f32, Y: mut &Vec<Vec<f32>>, ystep: mut &Vec<Vec<f32>>, gain: mut &Vec<Vec<f32>>, pij:  &Vec<Vec<f32>>
        let tco: f32 = optistep( n, dim, opticount, epsilon, &mut solY, &mut ystep, &mut gain, &pij, lastcost );
        //println!( "step {:?} cost {:?} sol {:?}", opticount, tco, solY );
        lastcost = tco;
        if opticount > maxopti {
            keeponopti = false;
        }
        opticount += 1;
    }
    //println!( "Y after opt: {:?}", solY );

    return solY;
}

/*---------------3D Data preparation and ...----------------------------------*/
fn dothelocomotion( obj_file: &str, nna:usize ) -> (Vec<f32>, usize) {
    //input a object file
    //FROM https://archaeologydataservice.ac.uk/catalogue/adsdata/arch-3369-1/dissemination/archaide_ceramics/obj/AM_72_DR591.obj
    //let obj_file = "/home/khk/Dokumente/scherbenmaschine/WS21-22/3ddate/Van_der_Werff_1_DR143.obj";
    //let obj_file = "/home/khk/Dokumente/scherbenmaschine/WS21-22/3ddate/Rhodian_Type_DR123.obj";

    let modelmy = tobj::load_obj( 
                    obj_file,
                    &tobj::LoadOptions{
                        triangulate: true,
                        single_index: true,
                        ..Default::default()
                    }, //checkout the options there is something for GPU usage
                  ).expect( "Failed to load OBJ from file!" );

    let momo = modelmy.0;

    let momomesh = &momo[0].mesh; //get the first and only mesh from the data

    //let i = momomesh.indices[0] as usize; //get index of WHAT - triangle first point 
    //let ii = momomesh.indices[1] as usize;
    //let facefucker = momomesh.face_arities.len();
    //println!( "Faces per mesh: {}", facefucker );
    
    //ALSO: indices array ist nicht identisch mit dem bloßen positzion array, die indices beziehen sich auf das vertex array, aber ich habe keine ahnung, wie die faces aussehen, jedenfalls kann ich die regelmaßigkeit nicht erkennen; since the index array len is a multiple of three (momomesh.indices.len()/3) I assume, that every three index bild a triangle 

    //divide the indices array into array of 3-vectors 
    //take a triangle and get all all neighbors, citeria 180 grad for a selected point

    //for krümmungsberechung musst du trjectories of points definieren und dies dann anhand der normalen auf den umgebenden Dreiecken untersuchen

    /*if facefucker == 0 { //computer your own faces
        // pos = [x, y, z]
        let pos_o = [
            
        ]; 
        let pos_oo = [
            momomesh.positions[ii * 3],
            momomesh.positions[ii * 3 + 1],
            momomesh.positions[ii * 3 + 2],
        ]; 
            // POS: [-2.527283, 774.58716, -181.54703], i: 0, ii: 1
        println!( "POS: {:?}, i: {}, POS: {:?}, ii: {}, indexarray/3: {:?}", pos_o, i, pos_oo, ii, momomesh.indices.len()/3 );
    }*/

    println!("Find extream vertices.");
    let mut pmaxz = [0.0,0.0,0.0]; //mutable arrays
    let mut pminz = [0.0,0.0,10000000.0];
    let mut pmaxx = [0.0,0.0,0.0];
    let mut pminx = [100000.0,0.0,0.0];
    let mut pmaxy = [0.0,0.0,0.0];
    let mut pminy = [0.0,1000000.0,0.0];
    
    for i in 0..momomesh.positions.len() / 3 {
        let p = [ momomesh.positions[i*3], momomesh.positions[i*3+1], momomesh.positions[i*3+2] ];
        if p[0] < pminx[0] {
            pminx[0] = p[0];
            pminx[1] = p[1];
            pminx[2] = p[2];
        }
        if p[1] < pminy[1] {
            pminy[0] = p[0];
            pminy[1] = p[1];
            pminy[2] = p[2];
        }
        if p[2] < pminz[2] {
            pminz[0] = p[0];
            pminz[1] = p[1];
            pminz[2] = p[2];
        }
        if p[0] > pmaxx[0] {
            pmaxx[0] = p[0];
            pmaxx[1] = p[1];
            pmaxx[2] = p[2];
        }
        if p[1] > pmaxy[1] {
            pmaxy[0] = p[0];
            pmaxy[1] = p[1];
            pmaxy[2] = p[2];
        }
        if p[2] > pmaxz[2] {
            pmaxz[0] = p[0];
            pmaxz[1] = p[1];
            pmaxz[2] = p[2];
        }
    }
    let mut basepointsofmodel = [pminx, pmaxx, pminy, pmaxy, pminz, pmaxz]; //subspace of modell pointspace
    //println!( "Basepoints: {:?}", basepointsofmodel );
    /*
    //get rotation axes and make it first point in array; we can do it like this bekause of rotsym of models
    if vec3len( &pminy ) - vec3len( &pmaxy ) != 0.0 {
        basepointsofmodel = [pminy, pmaxy, pminx, pmaxx, pminz, pmaxz];
    } else if vec3len( &pminz ) - vec3len( &pmaxz ) != 0.0 {
        basepointsofmodel = [pminz, pmaxz, pminx, pmaxx, pminy, pmaxy];
    }*/

    //build bounding box and take two seitenhalbierende along rot axes Y 
    let p1tm1 = [ (pmaxx[0].abs()-pminx[0].abs())/2.0, pmaxy[1], pmaxz[2]];
    let p2tm1 = [ pminx[0], pmaxy[1], pmaxz[2]];

    let p1t1 = [ (pmaxx[0].abs()-pminx[0].abs())/2.0, pmaxy[1], pmaxz[2]];
    let p2t1 = [ (pmaxx[0].abs()-pminx[0].abs())/2.0, pminy[1], pmaxz[2]];
    
    let p1t2 = [ pmaxx[0], pmaxy[1], (pmaxz[2].abs()-pminz[2].abs())/2.0];
    let p2t2 = [ pmaxx[0], pminy[1], (pmaxz[2].abs()-pminz[2].abs())/2.0];

    //println!( "Boundingbox Flächenhalbierende: p1t1 {:?} p2t1 {:?} und p1t2 {:?} p2t2 {:?}", p1t1, p2t1, p1t2, p2t2 );

    //build k-nearest on the model for the sampled pointes on tangent
    let tang0 = p32vec3( &p1tm1, &p2tm1 );
    let tang0l = vec3len( &tang0 );

    let tang1 = p32vec3( &p1t1, &p2t1 );
    let tang1l = vec3len( &tang1 );

    let tang2 = p32vec3( &p1t2, &p2t2 );
    let tang2l = vec3len( &tang2 );
    
    println!("Tangent len 0 {:?}, len 1 {:?}, len 2 {:?}", tang0l, tang1l, tang2l);

    //tangdesic: build subset of model with minimal distance to the 5 tangent
    //let mut d1 = Vec::new( );
    //let mut d2 = Vec::new( );
    let mut pointstodetermminpery = Vec::new( );
    for vtxi in 0..momomesh.indices.len( ) { // take all points of modell
        let i = momomesh.indices[ vtxi ] as usize; 
        let p = [ momomesh.positions[i], momomesh.positions[i+1], momomesh.positions[i+2], i as f32 ];
        pointstodetermminpery.push( p );
    }
    pointstodetermminpery.sort_by( |a, b| a[1].partial_cmp( &b[1] ).unwrap( ) );
    let mut oldyniveau = pointstodetermminpery[0][1];
    let mut yminimalpoints1 = Vec::new( );
    let mut yminimalpoints2 = Vec::new( );
    let mut currmindistd1 = 100000000000.0;
    let mut pointat1 = pointstodetermminpery[0];
    let mut currmindistd2 = 100000000000.0;
    let mut pointat2 = pointstodetermminpery[0];
    let mut selectedalongtangent1 = Vec::new( );
    let mut selectedalongtangent2 = Vec::new( );
    let mut niveausize = 0;
    for pp in 0..pointstodetermminpery.len() {
        if pointstodetermminpery[pp][1] != oldyniveau {
            if niveausize > 10 {
                yminimalpoints1.push(pointat1);
                selectedalongtangent1.push(pointat1[3] as usize);
                yminimalpoints2.push(pointat2);
                selectedalongtangent2.push(pointat2[3] as usize);
            } 
            currmindistd1 = 100000000000.0;
            currmindistd2 = 100000000000.0;
            niveausize = 0;
        }
        let temppoi = [pointstodetermminpery[pp][0], pointstodetermminpery[pp][1], pointstodetermminpery[pp][2]];
        let tempd1 = compdistp3vec3( &p2t1, &temppoi, &tang1, &tang1l );
        let tempd2 = compdistp3vec3( &p2t2, &temppoi, &tang2, &tang2l );
        if currmindistd1 > tempd1 {
            currmindistd1 = tempd1;
            pointat1 = pointstodetermminpery[pp];
        }
        if currmindistd2 > tempd2 {
            currmindistd2 = tempd2;
            pointat2 = pointstodetermminpery[pp];
        }
        oldyniveau = pointstodetermminpery[pp][1];
        niveausize += 1;
    }
    
    /*for vtxi in 0..momomesh.indices.len( ) { // take all points of modell
        let i = momomesh.indices[ vtxi ] as usize; 
        let p = [ momomesh.positions[i], momomesh.positions[i+1], momomesh.positions[i+2] ];
        //compare with gedesic1
        if( iteminarray( &yminimalpoints1, &p) ) {
            d1.push( compdistp3vec3( &p2t1, &p, &tang1, &tang1l ) );
            selectedalongtangent1.push(i);
        }
        if( iteminarray( &yminimalpoints2, &p) ) {
            d2.push( compdistp3vec3( &p2t2, &p, &tang2, &tang2l ) );
            selectedalongtangent2.push(i);
        }
    }
    println!("All distances to tangent 1 {:?}", d1.len() );*/
    //sort float vec
    /*let d1o = d1.clone();
    let d2o = d2.clone();
    d1.sort_by( |a, b| a.partial_cmp( b ).unwrap( ) ); //ist das aufsteigend oder absteigend???
    d2.sort_by( |a, b| a.partial_cmp( b ).unwrap( ) );
    
    //d1.iter_mut( ).for_each( |a, b| b - a );
    //let mut firstmaxindex = 0;
    let (absmaxindexd1, mut absmaxatd1) = getmaxdiffpos( &d1 );

    println!("theshold {:?}--{:?}", absmaxatd1, pmaxx[0]);
    if absmaxatd1 > pmaxx[0] {
        absmaxatd1 = (pmaxx[0]-(absmaxatd1 - pmaxx[0])).abs();
    }
    
    println!("D1 Max abs neg: {:?} with {:?}, total len: {:?}", absmaxindexd1, absmaxatd1, d1.len());
    let (absmaxindexd2, mut absmaxatd2) = getmaxdiffpos( &d2 );
    if absmaxatd2 > pmaxx[0] {
        absmaxatd2 = (pmaxx[0]-(absmaxatd2 - pmaxx[0])).abs();
    }
    println!("D2 Max abs neg: {:?} with {:?}, total len: {:?}", absmaxindexd2, absmaxatd2, d2.len());
    
    //select the d array with smallest absnaxindex
    /*let absmaxes = [absmaxatd1, absmaxatd2, absmaxatd3, absmaxatd4, absmaxatd5];
    let desall = [&d1, &d2, &d3, &d4, &d5];
    let mut so = absmaxes.clone();
    so.sort_by( |a, b| a.partial_cmp( b ).unwrap( ) );
    let mut des = Vec::new( );
    for i in 0..3 {
        for j in 0..absmaxes.len() {
            if absmaxes[j] == so[i] {
                println!("i: {:?} with {:?}, j: {:?} with {:?}", i, so[i], j, absmaxes[j]);
                des.push( desall[ i ] );
                break;
            }
        }
    }*/

    
    let mut cou = 0;
    let lim = 1000;//100000;
    for i in 0..d1o.len() {
        if lim <= cou {
            break;
        }
        //if d1o[i] < absmaxatd1*2.0 {
            selectedalongtangent1.push( i );
            cou += 1;
        //}
    }
    
    cou = 0;
    for i in 0..d2o.len() {
        if lim <= cou {
            break;
        }
        //if d2o[i] < absmaxatd2*2.0 {
            selectedalongtangent2.push( i );
            cou += 1;
        //}
    }
    */
    //println!("Selected indices related to tangent 1: {:?} of {:?}", selectedalongtangent1.len(), pointstodetermminpery.len() );
    //println!("Selected indices related to tangent 2: {:?} of {:?}", selectedalongtangent2.len(), pointstodetermminpery.len() );
    

    //println!( "Itterations for building triangles: {:?}", momomesh.indices.len()/3 );
    let mut triangles = Vec::new( );
    let mut pointtotriange = HashMap::new( );
    for vtxi in 0..momomesh.indices.len() / 3 { // build triangles
        let vv = vtxi*3;
        let i = momomesh.indices[ vv ] as usize; // p1 of triangle
        let ii = momomesh.indices[ vv+1 ] as usize; // p2 of triangle
        let iii = momomesh.indices[ vv+2 ] as usize; // p3 of triangle
        
        let p1 = [ momomesh.positions[i], momomesh.positions[i+1], momomesh.positions[i+2] ];
        let p1s = format!("{:.6}-{:.6}-{:.6}", momomesh.positions[i], momomesh.positions[i+1], momomesh.positions[i+2]);
        let p2 = [ momomesh.positions[ii], momomesh.positions[ii+1], momomesh.positions[ii+2] ];
        let p2s = format!("{:.6}-{:.6}-{:.6}", momomesh.positions[ii], momomesh.positions[ii+1], momomesh.positions[ii+2]);
        let p3 = [ momomesh.positions[iii], momomesh.positions[iii+1], momomesh.positions[iii+2] ];
        let p3s = format!("{:.6}-{:.6}-{:.6}", momomesh.positions[iii], momomesh.positions[iii+1], momomesh.positions[iii+2]);
        triangles.push( [ p1,p2,p3 ] );
        pointtotriange.entry( p1s ).or_insert( vec![] ).push( triangles.len()-1 );
        pointtotriange.entry( p2s ).or_insert( vec![] ).push( triangles.len()-1 );
        pointtotriange.entry( p3s ).or_insert( vec![] ).push( triangles.len()-1 );
    }
    //println!( "triangle 1: {:?}, triangle 2: {:?}, triangle 3: {:?}", triangles[0], triangles[1], triangles[2] );

    //nimm für jeden Punkt ein lattice an - das ist auch irgendwie UNGEIL, vielleicht mal nur eines nehmen
    let mut triangstars = Vec::new( );
    let mut starcp = Vec::new( );

    //first tangente is that with the henkel
    /*for vtxi in 0..selectedalongtangent1.len() {
        let i = momomesh.indices[ selectedalongtangent1[vtxi] ] as usize; //because - kann ich aus dem index eventuell schon ein lattice ableiten, ja kann sein, wenn im index drum rum mehr stellen den gleichen idex haben, das geht denke ich, jedenfalls muss man nicht über ein globes triangle gehen, das spart Zeit.
        let p = [ momomesh.positions[i], momomesh.positions[i+1], momomesh.positions[i+2] ];
        let ps = format!("{:.6}-{:.6}-{:.6}", momomesh.positions[i], momomesh.positions[i+1], momomesh.positions[i+2]);
        let mut lett = Vec::new();
        for r in 0..pointtotriange[ &ps ].len( ) {
            lett.push( triangles[ pointtotriange[ &ps ][ r ] ] );
        }
        triangstars.push( Tristar::new( p, lett ) );
        starcp.push( p );
    }*/
    //second tangent is that with the body
    for vtxi in 0..selectedalongtangent2.len() {
        let i = momomesh.indices[ selectedalongtangent2[vtxi] ] as usize; //because
        let p = [ momomesh.positions[i], momomesh.positions[i+1], momomesh.positions[i+2] ];
        let ps = format!("{:.6}-{:.6}-{:.6}", momomesh.positions[i], momomesh.positions[i+1], momomesh.positions[i+2]);
        let mut lett = Vec::new();
        for r in 0..pointtotriange[ &ps ].len( ) {
            lett.push( triangles[ pointtotriange[ &ps ][ r ] ] );
        }
        triangstars.push( Tristar::new( p, lett ) );
        starcp.push( p );
    }
    //sort triangle stars
    triangstars.sort_by( |a, b| a.poi[1].partial_cmp( &b.poi[1] ).unwrap( ) ); //SORTING DONE RIGHT?????? NEED TO sort by first index of first tupel partner

    //println!( "Trianglestars build: {:?}, len points {:?}, verhaeltnis ausgew sterne points {:?}", triangstars.len(), momomesh.indices.len(), (triangstars.len() as f32 / momomesh.indices.len() as f32)*1.0 );

    //build compare vector
    let mut compvec = Vec::new();
    let mut yisnow = triangstars[ 0 ].poi[ 1 ];
    let mut sumerg = 0.0;
    let mut countyrange = 0;
    for stari in 0..triangstars.len( ) {
        //println!("Lettice: {:?} with {:?} triangles, with curvature: {:?}", stari, triangstars[stari].len( ), angeldefect( &triangstars[stari], &starcp[stari] ) );
        sumerg += angeldefect( &triangstars[stari].sta, &triangstars[stari].poi );
        
        if yisnow != triangstars[stari].poi[1] {
            compvec.push( sumerg / countyrange as f32 );
            sumerg = 0.0;
            countyrange = 0;
        } 
        
        countyrange += 1;
        yisnow = triangstars[stari].poi[1];
    }

    /*for c in 0..compvec.len() {
        println!("{:?} -- {:?}", c, compvec[c] );
    }*/
    println!("comparation vector build");

    

    return (compvec, nna);
}


/*----------- INPUT OUTPUT of results ----------------------------------------*/
fn writejson( tocompare: &Vec<Vec<f32>>, namensind: &Vec<String>, nn:&str ) -> std::io::Result<()> {
    let mut bufferv = File::create(nn)?;
    let tw = (tocompare, namensind);
    let saverv = serde_json::to_writer_pretty(&bufferv, &tw).unwrap();
    //let mut buffern = File::create("nam.json")?;
    //let savenn = serde_json::to_writer(&buffern, &namensind).unwrap();
    Ok(saverv)

}

fn readjson(path: &Path) -> Result<( Vec<Vec<f32>>, Vec<String> ), Box<Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `User`.
    let u = serde_json::from_reader(reader)?;

    // Return the `User`.
    Ok(u)
}

/*---------------------running the main main main-----------------------------*/
fn main( ){
    //if needed
    let workinigdir = env::current_dir().unwrap();
    let wd = workinigdir.to_str().unwrap();
    println!("Me is the ricci revange by Dr. Ben running from {} ", workinigdir.display());
    

    //get data from channel
    let mut tocompare = Vec::new( );
    let mut namensind = Vec::new( );
    let fpstr = format!("{}{}", wd, "/vec.json");
    let fp = Path::new( &fpstr );
    if !fp.exists() {
        //let fp = fs::read_dir("/mnt/data/khk/scherbenbringen").unwrap();
        let fp = fs::read_dir("/home/khk/Dokumente/scherbenmaschine/WS21-22/3ddate/sm").unwrap();
        //let fp = fs::read_dir("/home/khk/Dokumente/scherbenmaschine/WS21-22/3ddate/simple").unwrap();
        let (tx, rx) = mpsc::channel();
        let mut hiha = Vec::new( );
     
        let cpus = num_cpus::get();
        println!("Num of CPU/cores {}", cpus );

        let mut spaned = 0;
        let mut tocomparena = Vec::new( );

        for p in fp {
            let entry = p.unwrap();
            let entry_path = entry.path();
            let file_name = entry_path.file_name().unwrap();
            let file_name_as_str = file_name.to_str().unwrap();
            let file_name_as_string = String::from(file_name_as_str);
        
            if file_name_as_string.contains( ".obj" ) {
                println!("Name: {}, {}", entry_path.display(), file_name_as_string );
                let txc = tx.clone();
                let cn = spaned.clone();
                tocomparena.push( file_name_as_string );
                let ha = thread::spawn( move || {
                    let past = entry_path.to_str().unwrap( );
                    txc.send( dothelocomotion( past, cn ) ).unwrap();
                    drop( txc );
                });
        
                hiha.push( ha );
                spaned += 1;
            }
        }
         
        //wait for the threads to complete
        for ha in hiha {
            ha.join().expect("Child thread panicked");
        }
        
        for _ in 0..spaned {
            let (cv, na) = rx.recv( ).unwrap( ); //give back name (filename) of related model
            println!("Vec to compare (len): {:?}", cv.len( ) ); //vec of vec floats
            tocompare.push( cv );
            namensind.push( tocomparena[ na ].clone() );
        }
        
        //serializeresults
        let wr =  writejson(&tocompare, &namensind, "vec.json");
    } else {
        //load result data of geometry computation from file
        let u = readjson(&fp).unwrap();
        println!( "Geom. data already present...");
        tocompare = u.0;
        namensind = u.1;
    }    


    //tSNE
    println!("tSNE doing with {:?} vectors of diff len.", tocompare.len( ) );

    let arrangedinlessd = tsne( &tocompare, 2, 20 ); //dim to embed to, number of neighbors, .... ; should give a drawable array

    writejson(&arrangedinlessd, &namensind, "resu.json");
    println!("Embedding {:?} Names {:?}", arrangedinlessd, namensind );
    //timing
    //save vec to folder
    println!( "End" );
}

