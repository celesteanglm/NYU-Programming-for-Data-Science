import time
"""
    N-body simulation.
    Change made: All Changes
    Time taken: 19.87
    R: 3.83
"""

def advance(BODIES, dt, iterations, bodies_keys, bodies_tbl):
    '''
        advance the system timestep by timestep
    '''
    # Reduce function call overhead by running the iterations within this function
    for _ in range(iterations):
        # Using dictionary instead of list for lookup
        seenit = dict(bodies_tbl)
        for body1 in bodies_keys:
            # Calling the 2 lines below within this inner loop to reduce runtime
            ([x1, y1, z1], v1, m1) = BODIES[body1]
            seenit[body1] = True
            for body2 in BODIES.keys():
                if not (seenit[body2]):
                    # Reduce function call overhead by moving below calculations here
                    # compute deltas
                    ([x2, y2, z2], v2, m2) = BODIES[body2]
                    (dx, dy, dz) = (x1-x2, y1-y2, z1-z2)
                    mag = dt * ((dx * dx + dy * dy + dz * dz) ** (-1.5))
                    # compute mag
                    m2_mag = m2 * mag
                    m1_mag = m1 * mag
                    # update vs
                    v1[0] -= dx * m2_mag
                    v1[1] -= dy * m2_mag
                    v1[2] -= dz * m2_mag
                    v2[0] += dx * m1_mag
                    v2[1] += dy * m1_mag
                    v2[2] += dz * m1_mag
            
        for body in bodies_keys:
            (r, [vx, vy, vz], m) = BODIES[body]
            # update rs
            r[0] += dt * vx
            r[1] += dt * vy
            r[2] += dt * vz
    
def report_energy(BODIES, bodies_keys, bodies_tbl, e=0.0):
    '''
        compute the energy and return it so that it can be printed
    '''
    # Using dictionary instead of list for lookup
    seenit = dict(bodies_tbl)
    for body1 in bodies_keys:
        # Calling the 2 lines below within this inner loop to reduce runtime
        ([x1, y1, z1], v1, m1) = BODIES[body1]
        seenit[body1] = True
        for body2 in bodies_keys:
            if not (seenit[body2]):
                ((x2, y2, z2), v2, m2) = BODIES[body2]
                (dx, dy, dz) = (x1-x2, y1-y2, z1-z2)
                # Reduce function call overhead by moving calculation here
                e -= (m1 * m2) / ((dx * dx + dy * dy + dz * dz) ** 0.5)
        
    for body in bodies_keys:
        (r, [vx, vy, vz], m) = BODIES[body]
        # Reduce function call overhead by moving calculation here
        e += m * (vx * vx + vy * vy + vz * vz) / 2.
        
    return e

def offset_momentum(BODIES, ref, bodies_keys, px=0.0, py=0.0, pz=0.0):
    '''
        ref is the body in the center of the system
        offset values from this reference
    '''
    for body in bodies_keys:
        (r, [vx, vy, vz], m) = BODIES[body]
        px -= vx * m
        py -= vy * m
        pz -= vz * m
        
    (r, v, m) = ref
    v[0] = px / m
    v[1] = py / m
    v[2] = pz / m


def nbody(loops, reference, iterations):
    '''
        nbody simulation
        loops - number of loops to run
        reference - body at center of system
        iterations - number of timesteps to advance
    '''
    # Using local instead of global variables
    PI = 3.14159265358979323
    SOLAR_MASS = 4 * PI * PI
    DAYS_PER_YEAR = 365.24

    BODIES = {
        'sun': ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], SOLAR_MASS),

        'jupiter': ([4.84143144246472090e+00,
                     -1.16032004402742839e+00,
                     -1.03622044471123109e-01],
                    [1.66007664274403694e-03 * DAYS_PER_YEAR,
                     7.69901118419740425e-03 * DAYS_PER_YEAR,
                     -6.90460016972063023e-05 * DAYS_PER_YEAR],
                    9.54791938424326609e-04 * SOLAR_MASS),

        'saturn': ([8.34336671824457987e+00,
                    4.12479856412430479e+00,
                    -4.03523417114321381e-01],
                   [-2.76742510726862411e-03 * DAYS_PER_YEAR,
                    4.99852801234917238e-03 * DAYS_PER_YEAR,
                    2.30417297573763929e-05 * DAYS_PER_YEAR],
                   2.85885980666130812e-04 * SOLAR_MASS),

        'uranus': ([1.28943695621391310e+01,
                    -1.51111514016986312e+01,
                    -2.23307578892655734e-01],
                   [2.96460137564761618e-03 * DAYS_PER_YEAR,
                    2.37847173959480950e-03 * DAYS_PER_YEAR,
                    -2.96589568540237556e-05 * DAYS_PER_YEAR],
                   4.36624404335156298e-05 * SOLAR_MASS),

        'neptune': ([1.53796971148509165e+01,
                     -2.59193146099879641e+01,
                     1.79258772950371181e-01],
                    [2.68067772490389322e-03 * DAYS_PER_YEAR,
                     1.62824170038242295e-03 * DAYS_PER_YEAR,
                     -9.51592254519715870e-05 * DAYS_PER_YEAR],
                    5.15138902046611451e-05 * SOLAR_MASS)}

    bodies_keys = BODIES.keys()
    bodies_tbl = dict([(body, False) for body in bodies_keys])
    # Set up global state
    offset_momentum(BODIES, BODIES[reference], bodies_keys)

    for _ in range(loops):
        report_energy(BODIES, bodies_keys, bodies_tbl)
        advance(BODIES, 0.01, iterations, bodies_keys, bodies_tbl)
        print(report_energy(BODIES, bodies_keys, bodies_tbl))

if __name__ == '__main__':
    start_time = time.time()
    nbody(100, 'sun', 20000)
    end_time = time.time()
    time_taken = end_time - start_time
    relative_speedup = 76.05/float(time_taken)
    print('Time taken: {0:.2f}'.format(time_taken))
    print('R: {0:.2f}'.format(relative_speedup))