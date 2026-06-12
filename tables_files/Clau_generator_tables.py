import mpmath as mp

mp.mp.dps = 80


# ============================================================
# Archivo de entrada y salida
# ============================================================

INPUT_FILE = "output_1.txt"
OUTPUT_FILE = "comparative_tables.txt"


# ============================================================
# Lectura del archivo por bloques
# ============================================================

def read_blocks(filename):
    blocks = []
    header = None
    data = []

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                if header is not None:
                    blocks.append((header, data))
                header = line
                data = []
            else:
                data.append(line.split())

    if header is not None:
        blocks.append((header, data))

    if len(blocks) != 5:
        raise ValueError(
            f"Se esperaban 5 bloques de datos, pero se encontraron {len(blocks)}."
        )

    return blocks


def mp_column(data, col):
    return [mp.mpf(row[col]) for row in data]


def parse_file(filename):
    blocks = read_blocks(filename)

    # Bloque 1: SciPy/Numpy, 16 dígitos: t, n(t)
    t = mp_column(blocks[0][1], 0)
    n_numpy = mp_column(blocks[0][1], 1)

    # Bloque 2: Runge-Kutta, 32 dígitos: n(t)
    n_rk = mp_column(blocks[1][1], 0)

    # Bloque 3: Runge-Kutta, 32 dígitos: C(t)
    C_rk = mp_column(blocks[2][1], 0)

    # Bloque 4: método propio con mpmath, 32 dígitos: t, n(t)
    t_mpmath = mp_column(blocks[3][1], 0)
    n_mpmath = mp_column(blocks[3][1], 1)

    # Bloque 5: método propio con SciPy/Numpy, 16 dígitos: t, C(t)
    t_C = mp_column(blocks[4][1], 0)
    C_numpy = mp_column(blocks[4][1], 1)

    lengths = [
        len(t),
        len(n_numpy),
        len(n_rk),
        len(C_rk),
        len(t_mpmath),
        len(n_mpmath),
        len(t_C),
        len(C_numpy),
    ]

    if len(set(lengths)) != 1:
        raise ValueError(f"Los bloques no tienen la misma longitud: {lengths}")

    if t != t_mpmath:
        raise ValueError("Los tiempos de n(t) con mpmath no coinciden.")

    if t != t_C:
        raise ValueError("Los tiempos de C(t) no coinciden.")

    return t, n_numpy, n_rk, C_rk, n_mpmath, C_numpy


# ============================================================
# Formato numérico
# ============================================================

def fmt_time(x):
    return str(int(x))


def fmt_number(x, digits):
    return mp.nstr(x, digits)


def ape(reference, approx):
    if reference == 0:
        return mp.nan
    return abs((approx - reference) / reference) * 100


def fmt_ape(x, digits=3):
    if mp.isnan(x):
        return "nan"
    return mp.nstr(x, digits)


# ============================================================
# Escritura de tablas
# ============================================================

def write_table(
    f,
    title,
    t_values,
    rk_values,
    reference_values,
    rk_digits,
    reference_digits,
    only_even=False,
    ape_digits=3,
):
    f.write(title + "\n")
    f.write("-" * len(title) + "\n")

    f.write(
        f"{'t':>4} "
        f"{'Runge-Kutta':>42} "
        f"{'Metodo propio / referencia':>42} "
        f"{'APE (%)':>12}\n"
    )

    for t, rk, ref in zip(t_values, rk_values, reference_values):
        if only_even and int(t) % 2 != 0:
            continue

        error = ape(ref, rk)

        f.write(
            f"{fmt_time(t):>4} "
            f"{fmt_number(rk, rk_digits):>42} "
            f"{fmt_number(ref, reference_digits):>42} "
            f"{fmt_ape(error, ape_digits):>12}\n"
        )

    f.write("\n\n")


def create_tables(input_file, output_file):
    t, n_numpy, n_rk, C_rk, n_mpmath, C_numpy = parse_file(input_file)

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("=" * 90 + "\n")
        f.write(f"TABLAS COMPARATIVAS PARA: {input_file}\n")
        f.write("=" * 90 + "\n\n")

        # ====================================================
        # Tablas con todos los tiempos
        # ====================================================

        write_table(
            f,
            "Tabla 1. n(t): Runge-Kutta vs SciPy/Numpy",
            t,
            n_rk,
            n_numpy,
            rk_digits=32,
            reference_digits=16,
            only_even=False,
            ape_digits=3,
        )

        write_table(
            f,
            "Tabla 2. n(t): Runge-Kutta vs metodo propio con mpmath",
            t,
            n_rk,
            n_mpmath,
            rk_digits=32,
            reference_digits=32,
            only_even=False,
            ape_digits=4,
        )

        write_table(
            f,
            "Tabla 3. C(t): Runge-Kutta vs metodo propio con SciPy/Numpy",
            t,
            C_rk,
            C_numpy,
            rk_digits=32,
            reference_digits=16,
            only_even=False,
            ape_digits=3,
        )

        # ====================================================
        # Tablas solo con tiempos pares
        # ====================================================

        write_table(
            f,
            "Tabla 4. n(t): Runge-Kutta vs SciPy/Numpy, tiempos pares",
            t,
            n_rk,
            n_numpy,
            rk_digits=32,
            reference_digits=16,
            only_even=True,
            ape_digits=3,
        )

        write_table(
            f,
            "Tabla 5. n(t): Runge-Kutta vs metodo propio con mpmath, tiempos pares",
            t,
            n_rk,
            n_mpmath,
            rk_digits=32,
            reference_digits=32,
            only_even=True,
            ape_digits=4,
        )

        write_table(
            f,
            "Tabla 6. C(t): Runge-Kutta vs metodo propio con SciPy/Numpy, tiempos pares",
            t,
            C_rk,
            C_numpy,
            rk_digits=32,
            reference_digits=16,
            only_even=True,
            ape_digits=3,
        )


# ============================================================
# Ejecución
# ============================================================

if __name__ == "__main__":
    create_tables(INPUT_FILE, OUTPUT_FILE)

    print(f"Archivo creado correctamente: {OUTPUT_FILE}")
