const N = 64;
const iter = 4;
const SCALE = 4;

function IX(x, y)
{
    return x + y * N;
}

function set_bnd(b, x)
{
    for(let i = 1; i < N - 1; i++) {
        x[IX(i, 0  )] = b == 2 ? -x[IX(i, 1  )] : x[IX(i, 1  )];
        x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
    }
    for(let j = 1; j < N - 1; j++) {
        x[IX(0  , j)] = b == 1 ? -x[IX(1  , j)] : x[IX(1  , j)];
        x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
    }

    x[IX(0, 0)]       = 0.33 * (x[IX(1, 0)] + x[IX(0, 1)] + x[IX(0, 0)]);
    x[IX(0, N-1)]     = 0.33 * (x[IX(1, N-1)] + x[IX(0, N-2)] + x[IX(0, N-1)]);
    x[IX(N-1, 0)]     = 0.33 * (x[IX(N-2, 0)] + x[IX(N-1, 1)] + x[IX(N-1, 0)]);
    x[IX(N-1, N-1)]   = 0.33 * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)] + x[IX(N-1, N-1)]);
}


function lin_solve(b, x, x0, a, c)
{
    const cRecip = 1 / c;
    for (let k = 0; k < iter; k++) {
        for (let j = 1; j < N - 1; j++) {
            for (let i = 1; i < N - 1; i++) {
                x[IX(i, j)] =
                    (x0[IX(i, j)]
                        + a*(x[IX(i+1, j  )]
                            +x[IX(i-1, j  )]
                            +x[IX(i  , j+1)]
                            +x[IX(i  , j-1)]
                        )) * cRecip;
            }
        }
        set_bnd(b, x);
    }
}


function diffuse (b, x, x0, diff, dt)
{
    const a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
}

function project(velocX, velocY, p, div)
{
    for (let j = 1; j < N - 1; j++) {
        for (let i = 1; i < N - 1; i++) {
            div[IX(i, j)] = -0.5*(
                velocX[IX(i+1, j    )]
                -velocX[IX(i-1, j    )]
                +velocY[IX(i  , j+1  )]
                -velocY[IX(i  , j-1  )]
            )/N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div); 
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);

    for (let j = 1; j < N - 1; j++) {
        for (let i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] -= 0.5 * (  p[IX(i+1, j)]
                -p[IX(i-1, j)]) * N;
            velocY[IX(i, j)] -= 0.5 * (  p[IX(i, j+1)]
                -p[IX(i, j-1)]) * N;
        }
    }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}


function advect(b, d, d0, velocX, velocY, dt)
{
    let i0, i1, j0, j1;

    let dtx = dt * (N - 2);
    let dty = dt * (N - 2);

    let s0, s1, t0, t1;
    let tmp1, tmp2, x, y;

    let Nfloat = N - 2;
    let ifloat, jfloat;

    for(let j = 1; j < N - 1; j++) { 
        for(let i = 1; i < N - 1; i++) {
            tmp1 = dtx * velocX[IX(i, j)];
            tmp2 = dty * velocY[IX(i, j)];
            x    = i - tmp1; 
            y    = j - tmp2;

            if(x < 0.5) x = 0.5; 
            if(x > Nfloat + 0.5) x = Nfloat + 0.5; 
            i0 = Math.floor(x); 
            i1 = i0 + 1;
            if(y < 0.5) y = 0.5; 
            if(y > Nfloat + 0.5) y = Nfloat + 0.5; 
            j0 = Math.floor(y);
            j1 = j0 + 1; 

            s1 = x - i0; 
            s0 = 1.0 - s1; 
            t1 = y - j0; 
            t0 = 1.0 - t1;

            d[IX(i, j)] = 
                s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }
    set_bnd(b, d);
}


class Fluid
{
    size;
    dt;
    diff;
    visc;

    s; density;
    Vx; Vy;
    Vx0; Vy0; 

    diff(val)
    {
        this.diffuse = val;
    }

    visc(val)
    {
        this.viscosity = val;
    }

    step(val)
    {
        this.dt = val;
    }

    constructor(diffusion, viscosity, dt)
    {
        this.size = N;
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;

        this.s = new Array(N*N).fill(0);
        this.density = new Array(N*N).fill(0);

        this.Vx = new Array(N*N).fill(0);
        this.Vy = new Array(N*N).fill(0);

        this.Vx0 = new Array(N*N).fill(0);
        this.Vy0 = new Array(N*N).fill(0);
    }

    cubeStep ()
    {
        let visc     = this.visc;
        let diff     = this.diff;
        let dt       = this.dt;

        let Vx      = this.Vx;
        let Vy      = this.Vy;

        let Vx0     = this.Vx0;
        let Vy0     = this.Vy0;

        let s       = this.s;
        let density = this.density;

        diffuse(1, Vx0, Vx, visc, dt);
        diffuse(2, Vy0, Vy, visc, dt);

        project(Vx0, Vy0, Vx, Vy);

        advect(1, Vx, Vx0, Vx0, Vy0, dt);
        advect(2, Vy, Vy0, Vx0, Vy0, dt);

        project(Vx, Vy, Vx0, Vy0);

        diffuse(0, s, density, diff, dt);
        advect(0, density, s, Vx, Vy, dt);
    }

    addDensity(x, y, amount)
    {
        let index = IX(x, y);
        this.density[index] += amount;
    }   

    addVelocity(x, y, amountX, amountY)
    {
        let index = IX(x, y);
        this.Vx[index] += amountX;
        this.Vy[index] += amountY;
    }

    renderD ()
    {
        translate(-width/2, -height/2);

        for (let i = 0; i < N; i++)
        {
            for (let j = 0; j < N; j++)
            {
                let x = i * SCALE;
                let y = j * SCALE;
                let d = this.density[IX(i, j)];
                fill(d);
                noStroke();
                rect(x, y, SCALE, SCALE);
            }
        }
    }

    fadeD ()
    {
        for (let i = 0; i < this.density.length; i++)
        {
            const d = this.density[i];
            this.density[i] = constrain(d-0.6, 0, 255);
        }
    }
}

let fluid;
let diff, visc, step;

function setup()
{
    createCanvas(N * SCALE, N * SCALE, WEBGL);
    diff = document.getElementById("diff");
    visc = document.getElementById("visc");
    step = document.getElementById("step");
    fluid = new Fluid(0, 0, 0);
}

function mouseDragged()
{
    let x = Math.floor(mouseX / SCALE);
    let y = Math.floor(mouseY / SCALE);
    if (x >= 0 && x < N && y >= 0 && y < N) {
        fluid.addDensity(x, y, 500);
        let amountX = (mouseX - pmouseX) * 2;
        let amountY = (mouseY - pmouseY) * 2;
        fluid.addVelocity(x, y, amountX, amountY);
    }

    diff.addEventListener("change", (e) => {
        fluid.diff(e.target.value);
    });

    visc.addEventListener("change", (e) => {
        fluid.visc(e.target.value);
    });

    step.addEventListener("change", (e) => {
        fluid.step(e.target.value);
    });
}

function draw()
{
    background(0);
    fluid.cubeStep();
    fluid.renderD();
    fluid.fadeD();
}
