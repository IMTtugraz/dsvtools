#!/usr/bin/python3

gx_unit=1e-3
dt_unit=1e-6
gamma=42.58e6

def normalizer(fov, gx_unit=gx_unit, gamma=gamma, dt_unit=dt_unit):
    #gx_unit * dt_unit * gamma * 2 * np.pi / (samples * 2*np.pi/fov) * samples
    return gx_unit * dt_unit * gamma * fov

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("fov", type=float, help="FOV (m)")
    parser.add_argument("--gx_unit", type=float, help="Gradient unit (T/m)", default=gx_unit)
    parser.add_argument("--dt_unit", type=float, help="ADC raster time (s)", default=dt_unit)
    parser.add_argument("--gamma", type=float, help="Gyromagnetic ratio", default=gamma)
    args = parser.parse_args()

    print(normalizer(args.fov, gx_unit=args.gx_unit, gamma=args.gamma, dt_unit=args.dt_unit))
