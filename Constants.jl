module Constants
    # Does not export anything - rely on namespace for use

    set_zero_subnormals(true)

    # Histogram and problem parameters
    global const num_t = convert(Int64, 1e5)
    global const num_particles = convert(Int64, 1e3)

    # Geometry definitions (assume minimum values are at origin)
    global const x_len = 1e-1  # cm
    global const y_len = 1e-1  # cm
    global const z_len = 1e-1  # cm
    global const vol = x_len * y_len * z_len  # cm^3

    # Point source position (at midpoint)
    global const src_x = x_len / 2.0  # cm
    global const src_y = y_len / 2.0  # cm
    global const src_z = z_len / 2.0  # cm

    # Initial conditions
    global const init_intensity = 1e0  # erg/cm^2-s
    global const init_temp = 1e0  # eV

    # Iteration condition
    global const num_say = convert(Int64, num_t / 1e2)
    global const weight_cutoff = 1e-10

    # Physics constants
    global const sol = 29979245800.0  # cm/s
    global const arad = 137.2017050419133  # erg/cm^3-eV^4
    global const sb_const = 5.6704e-5 * 11604.525^4  # erg/cm^2-s-eV^4

    # Problem physics parameters
    global const t_max = 1e-11  # s
    global const t_init = 0.0  # s
    global const delta_t = (t_max - t_init) / num_t  # s
    global const chord_1 = 1e-3 / sol  # s
    global const chord_2 = 5e-3 / sol  # s
    global const dens_1 = 1.0  # g/cm^3
    global const dens_2 = 1.0  # g/cm^3
    global const opacity_1 = 1.0  # cm^-1
    global const opacity_2 = 5.0  # cm^-1
    global const spec_heat_1 = 1.0  # erg/g-eV
    global const spec_heat_2 = 1.0  # erg/g-eV

    # Problem volumes
    global const volfrac_1 = chord_1 / (chord_1 + chord_2)
    global const volfrac_2 = 1.0 - volfrac_1
    global const vol_1 = volfrac_1 * vol  # cm^3
    global const vol_2 = volfrac_2 * vol  # cm^3
end
