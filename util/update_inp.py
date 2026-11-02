#!/usr/bin/env python3
import sys
from datetime import datetime

# Cutoff dates for historical input changes
D_2021_07_15 = datetime(2021, 7, 15)
D_2021_08_10 = datetime(2021, 8, 10)
D_2023_01_05 = datetime(2023, 1, 5)
D_2023_01_07 = datetime(2023, 1, 7)
D_2023_01_17 = datetime(2023, 1, 17)
D_2023_02_01 = datetime(2023, 2, 1)
D_2023_02_04 = datetime(2023, 2, 4)
D_2025_10_17 = datetime(2025, 10, 17)
D_2026_05_27 = datetime(2026, 5, 27)
D_2026_06_12 = datetime(2026, 6, 12)
D_2026_07_30 = datetime(2026, 7, 30)
D_2026_08_02 = datetime(2026, 8, 2)

def parse_date(date_str):
    for fmt in ('%Y-%m-%d', '%Y/%m/%d', '%Y%m%d'):
        try:
            return datetime.strptime(date_str.strip(), fmt)
        except ValueError:
            pass
    raise ValueError(f"Could not parse date: '{date_str}'. Please use YYYY-MM-DD format.")

def get_clean_lines(filename):
    clean_lines = []
    preceding_comments = []
    with open(filename, 'r') as f:
        for line in f:
            line_str = line.rstrip('\n')
            line_strip = line_str.strip()
            if not line_strip:
                preceding_comments.append((True, "")) # blank line
                continue
            if line_strip.startswith(';'):
                preceding_comments.append((False, line_str)) # comment line
                continue
            
            # This is a parameter line
            clean_content = line_strip
            trailing_comment = ""
            if ';' in line_strip:
                parts = line_strip.split(';', 1)
                clean_content = parts[0].strip()
                trailing_comment = ';' + parts[1]
                
            clean_lines.append({
                'content': clean_content,
                'trailing': trailing_comment,
                'preceding': preceding_comments
            })
            preceding_comments = []
    return clean_lines, preceding_comments

def parse_inp(clean_lines, date):
    ptr = 0
    comments_map = {}

    def next_line(key=None):
        nonlocal ptr
        if ptr >= len(clean_lines):
            raise ValueError(f"Unexpected end of input file at parameter line {ptr}")
        line_info = clean_lines[ptr]
        ptr += 1
        if key is not None:
            comments_map[key] = (line_info['preceding'], line_info['trailing'])
        return line_info['content']

    # Helper to get tokens of next line
    def get_tokens(key=None):
        line = next_line(key)
        return line.replace(',', ' ').split()

    # 1. Mesh parameters
    tokens = get_tokens('nex_nez')
    nex, nez = int(tokens[0]), int(tokens[1])
    
    tokens = get_tokens('x0_z0')
    x0, z0 = float(tokens[0]), float(tokens[1])
    
    tokens = get_tokens('rxbo_rzbo')
    rxbo, rzbo = float(tokens[0]), float(tokens[1])
    
    tokens = get_tokens('ircoord')
    ircoord = int(tokens[0])
    coordfile = tokens[1] if len(tokens) > 1 else ""

    tokens = get_tokens('nzonx')
    nzonx = int(tokens[0])
    zones_x = []
    for idx in range(nzonx):
        tk = get_tokens(f'zones_x_{idx}')
        zones_x.append((int(tk[0]), float(tk[1])))

    tokens = get_tokens('nzony')
    nzony = int(tokens[0])
    zones_y = []
    for idx in range(nzony):
        tk = get_tokens(f'zones_y_{idx}')
        zones_y.append((int(tk[0]), float(tk[1])))

    # 2. Tracer/Marker Section (old format only)
    if date < D_2023_01_17:
        if date < D_2021_07_15:
            # iint_marker, iint_tracer
            get_tokens('tracer_junk1')
            nzone_marker = int(get_tokens('tracer_junk2')[0])
            for idx in range(nzone_marker):
                get_tokens(f'tracer_zone_{idx}')
            # nzone_tracer, dt_outtracer
            get_tokens('tracer_junk3')
            for idx in range(nzone_marker):
                get_tokens(f'tracer_zone2_{idx}')
        else:
            # iint_marker, i_junk
            get_tokens('tracer_junk1')
            nzone_marker = int(get_tokens('tracer_junk2')[0])
            for idx in range(nzone_marker):
                get_tokens(f'tracer_zone_{idx}')
            # n_junk, d_junk
            get_tokens('tracer_junk3')

    # 3. Mechanical boundary conditions
    tokens = get_tokens('nystressbc_nydrsides')
    nystressbc, nydrsides = int(tokens[0]), int(tokens[1])
    
    tokens = get_tokens('nofbc')
    nofbc = int(tokens[0])
    bcs = []
    for idx in range(nofbc):
        tk = get_tokens(f'bc_{idx}')
        nofside = int(tk[0])
        nbc1 = int(tk[1])
        nbc2 = int(tk[2])
        nbc = int(tk[3])
        vals = list(map(float, tk[4:13]))
        bcs.append((nofside, nbc1, nbc2, nbc) + tuple(vals))

    tokens = get_tokens('nyhydro_pisos_iphsub_drosub_damp_vis')
    nyhydro = int(tokens[0])
    pisos = float(tokens[1])
    iphsub = int(tokens[2])
    drosub = float(tokens[3])
    damp_vis = float(tokens[4])

    tokens = get_tokens('g')
    g = float(tokens[0])

    # 4. Thermal
    tokens = get_tokens('i_prestress')
    i_prestress = int(tokens[0])

    if date >= D_2025_10_17:
        tokens = get_tokens('extra_pres')
        extra_pres = float(tokens[0])
    else:
        extra_pres = 0.0e9

    tokens = get_tokens('itherm')
    itherm = int(tokens[0])
    
    tokens = get_tokens('istress_therm')
    istress_therm = int(tokens[0])
    
    tokens = get_tokens('ishearh')
    ishearh = int(tokens[0])
    
    tokens = get_tokens('t_top')
    t_top = float(tokens[0])
    
    tokens = get_tokens('t_bot')
    t_bot = float(tokens[0])

    tokens = get_tokens('hs_hr')
    hs, hr = float(tokens[0]), float(tokens[1])
    
    tokens = get_tokens('itemp_bc_bot_bc')
    itemp_bc, bot_bc = int(tokens[0]), float(tokens[1])

    tokens = get_tokens('irtemp')
    irtemp = int(tokens[0])
    
    tokens = get_tokens('tempfile')
    tempfile = tokens[0]

    if date < D_2023_01_17:
        # time_scale
        get_tokens('time_scale_junk')

    tokens = get_tokens('nzone_age')
    nzone_age = int(tokens[0])
    zones_age = []
    for idx in range(nzone_age):
        first_line = get_tokens(f'zone_age_{idx}_line1')
        if len(first_line) >= 10:  # Old 1-line format
            if date < D_2023_01_05:
                ictherm = 1
                age_1 = float(first_line[0])
                tp1 = 0.0
                tp2 = 0.0
                hcs = list(map(float, first_line[1:5]))
                iphs = list(map(int, first_line[5:10]))
                ixtb1 = int(first_line[10])
                ixtb2 = int(first_line[11])
                nph_layer = 5
            else:
                ictherm = int(first_line[0])
                age_1 = float(first_line[1])
                tp1 = float(first_line[2])
                tp2 = float(first_line[3])
                hcs = list(map(float, first_line[4:8]))
                iphs = list(map(int, first_line[8:13]))
                ixtb1 = int(first_line[13])
                ixtb2 = int(first_line[14])
                nph_layer = 5
            zones_age.append((ictherm, age_1, tp1, tp2, nph_layer, hcs, iphs, ixtb1, ixtb2))
        else:  # 3-line format
            if len(first_line) == 7:  # nph_layer on line 1
                ictherm = int(first_line[0])
                age_1 = float(first_line[1])
                tp1 = float(first_line[2])
                tp2 = float(first_line[3])
                nph_layer = int(first_line[4])
                ixtb1 = int(first_line[5])
                ixtb2 = int(first_line[6])

                line2 = get_tokens(f'zone_age_{idx}_line2')
                hcs = list(map(float, line2))

                line3 = get_tokens(f'zone_age_{idx}_line3')
                iphs = list(map(int, line3))
            else:  # nph_layer prepended to line 2 (newest)
                ictherm = int(first_line[0])
                age_1 = float(first_line[1])
                tp1 = float(first_line[2])
                tp2 = float(first_line[3])
                ixtb1 = int(first_line[4])
                ixtb2 = int(first_line[5])

                line2 = get_tokens(f'zone_age_{idx}_line2')
                nph_layer = int(line2[0])
                hcs = list(map(float, line2[1:]))

                line3 = get_tokens(f'zone_age_{idx}_line3')
                iphs = list(map(int, line3))
            
            # Pad hcs and iphs to 5 elements for uniform representation
            hcs = hcs + [0.0] * (5 - len(hcs))
            iphs = iphs + [1] * (5 - len(iphs))
            zones_age.append((ictherm, age_1, tp1, tp2, nph_layer, hcs, iphs, ixtb1, ixtb2))

    # 5. Rheology
    tokens = get_tokens('nphase')
    nphase = int(tokens[0])
    phases = []
    for idx in range(nphase):
        tk = get_tokens(f'phase_{idx}')
        irheol = int(tk[0])
        visc = float(tk[1])
        den = float(tk[2])
        alfa = float(tk[3])
        beta = float(tk[4])
        pln = float(tk[5])
        acoef = float(tk[6])
        eactiv = float(tk[7])
        
        # Auto-detect vactiv based on parameter count
        if len(tk) >= 25:
            vactiv = float(tk[8])
            rest = list(map(float, tk[9:25]))
        else:
            vactiv = 0.0
            rest = list(map(float, tk[8:24]))

        phases.append((irheol, visc, den, alfa, beta, pln, acoef, eactiv, vactiv) + tuple(rest))

    # 6. Phase layers / horizontal layers (only if date < 2023-01-17)
    if date < D_2023_01_17:
        mphase = int(get_tokens('mphase_junk')[0])
        nphasl = int(get_tokens('nphasl_junk')[0])
        for idx in range(nphasl):
            get_tokens(f'phasl_{idx}_junk')

    # 7. Initial phase distribution
    tokens = get_tokens('irphase')
    irphase = int(tokens[0])
    phasefile = get_tokens('phasefile')[0]

    # 8. Inhomogeneities
    tokens = get_tokens('inhom')
    inhom = int(tokens[0])
    inhoms = []
    for idx in range(inhom):
        tk = get_tokens(f'inhom_{idx}')
        ix1_i = int(tk[0])
        ix2_i = int(tk[1])
        iy1_i = int(tk[2])
        iy2_i = int(tk[3])
        inphase_i = int(tk[4])
        igeom_i = int(tk[5])
        amp = float(tk[6])
        inhoms.append((ix1_i, ix2_i, iy1_i, iy2_i, inphase_i, igeom_i, amp))

    # 9. Limits & Healing
    tokens = get_tokens('ten_off')
    ten_off = float(tokens[0])
    
    tokens = get_tokens('tau_heal')
    tau_heal = float(tokens[0])
    
    v_lims = get_tokens('v_lims')
    v_min = float(v_lims[0])
    v_max = float(v_lims[1])
    ivis_shape = int(v_lims[2])
    efoldc = float(v_lims[3])

    # Melting parameters
    if date < D_2021_08_10:
        line1 = get_tokens('magma_line1')
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        prod_magma = float(line1[2])
        nelem_dike = 1
        prod_magma2 = 0.0
        rho_magma = 2700.0

        line2 = get_tokens('magma_line2')
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[3])
        ratio_crust_mzone = 0.1

        line3 = get_tokens('magma_line3')
        latent_heat_magma = 4.2e5
        lambda_freeze = float(line3[0])
        lambda_freeze_tdep = float(line3[1])

        line4 = get_tokens('magma_line4')
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    elif date < D_2023_02_01:
        line1 = get_tokens('magma_line1')
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        prod_magma = float(line1[2])
        nelem_dike = 1
        prod_magma2 = 0.0
        rho_magma = 2700.0

        line2 = get_tokens('magma_line2')
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[3])
        ratio_crust_mzone = 0.1

        line3 = get_tokens('magma_line3')
        latent_heat_magma = 4.2e5
        lambda_freeze = float(line3[0])
        lambda_freeze_tdep = float(line3[1])

        line4 = get_tokens('magma_line4')
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    elif date < D_2023_02_04:
        line1 = get_tokens('magma_line1')
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        prod_magma = float(line1[2])
        rho_magma = float(line1[3])
        nelem_dike = 1
        prod_magma2 = 0.0

        line2 = get_tokens('magma_line2')
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[3])
        ratio_crust_mzone = 0.1

        line3 = get_tokens('magma_line3')
        latent_heat_magma = float(line3[0])
        lambda_freeze = float(line3[1])
        lambda_freeze_tdep = float(line3[2])

        line4 = get_tokens('magma_line4')
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    elif date < D_2026_05_27:
        line1 = get_tokens('magma_line1')
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        prod_magma = float(line1[2])
        rho_magma = float(line1[3])
        nelem_dike = 1
        prod_magma2 = 0.0

        line2 = get_tokens('magma_line2')
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[2])
        ratio_crust_mzone = 0.1

        line3 = get_tokens('magma_line3')
        latent_heat_magma = float(line3[0])
        lambda_freeze = float(line3[1])
        lambda_freeze_tdep = float(line3[2])

        line4 = get_tokens('magma_line4')
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    elif date < D_2026_07_30:
        line1 = get_tokens('magma_line1')
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        nelem_dike = int(line1[2])
        prod_magma = float(line1[3])
        rho_magma = float(line1[4])
        prod_magma2 = 0.0

        line2 = get_tokens('magma_line2')
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[2])
        ratio_crust_mzone = 0.1

        line3 = get_tokens('magma_line3')
        latent_heat_magma = float(line3[0])
        lambda_freeze = float(line3[1])
        lambda_freeze_tdep = float(line3[2])

        line4 = get_tokens('magma_line4')
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    elif date < D_2026_08_02:
        line1 = get_tokens('magma_line1')
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        nelem_dike = int(line1[2])
        prod_magma = float(line1[3])
        prod_magma2 = float(line1[4])
        rho_magma = float(line1[5])

        line2 = get_tokens('magma_line2')
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[2])
        ratio_crust_mzone = 0.1

        line3 = get_tokens('magma_line3')
        latent_heat_magma = float(line3[0])
        lambda_freeze = float(line3[1])
        lambda_freeze_tdep = float(line3[2])

        line4 = get_tokens('magma_line4')
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    else:
        line1 = get_tokens('magma_line1')
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        nelem_dike = int(line1[2])
        prod_magma = float(line1[3])
        prod_magma2 = float(line1[4])
        rho_magma = float(line1[5])

        line2 = get_tokens('magma_line2')
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[2])
        ratio_crust_mzone = float(line2[3])

        line3 = get_tokens('magma_line3')
        latent_heat_magma = float(line3[0])
        lambda_freeze = float(line3[1])
        lambda_freeze_tdep = float(line3[2])

        line4 = get_tokens('magma_line4')
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    # 10. Remeshing
    rem_line = get_tokens('rem_line')
    ny_rem = int(rem_line[0])
    mode_rem = int(rem_line[1])
    ntest_rem = int(rem_line[2])
    angle_rem = float(rem_line[3])

    dx_rem = float(get_tokens('dx_rem')[0])
    
    tokens = get_tokens('topo_kappa_fac_kappa')
    topo_kappa, fac_kappa = float(tokens[0]), float(tokens[1])

    # 11. Process control
    idt_scale = int(get_tokens('idt_scale')[0])
    
    tokens = get_tokens('dt_scale_tolerance')
    dt_scale, tolerance = float(tokens[0]), float(tokens[1])
    
    tokens = get_tokens('i_rey_xReyn')
    i_rey, xReyn = int(tokens[0]), float(tokens[1])

    if_rmasses = int(get_tokens('if_rmasses')[0])
    if_imasses = int(get_tokens('if_imasses')[0])
    if_visc = int(get_tokens('if_visc')[0])
    if_avgsr = int(get_tokens('if_avgsr')[0])

    tokens = get_tokens('frac_elastic_frac_maxwell')
    frac_elastic, frac_maxwell = float(tokens[0]), float(tokens[1])
    
    tokens = get_tokens('movegrid_ndim')
    movegrid, ndim = int(tokens[0]), int(tokens[1])
    
    tokens = get_tokens('damping_mix1_mix2')
    damping, mix1, mix2 = float(tokens[0]), int(tokens[1]), int(tokens[2])

    # 12. Output
    time_max = float(get_tokens('time_max')[0])
    dtout_screen = float(get_tokens('dtout_screen')[0])
    dtout_file = float(get_tokens('dtout_file')[0])

    if date < D_2023_01_17:
        nprofil = int(get_tokens('nprofil_junk')[0])
        for idx in range(nprofil):
            get_tokens(f'profil_{idx}_junk')

    io_flags = list(map(int, get_tokens('io_flags')))
    lastout = int(get_tokens('lastout')[0])
    dtsave_file = float(get_tokens('dtsave_file')[0])
    lastsave = int(get_tokens('lastsave')[0])

    # Phase mapping for updates before 2026-06-12:
    # koceandry was 1, ksed1 was 10, ksed2 was 11.
    # Now koceandry is 10, ksed1 is 1, ksed2 is 11.
    if date < D_2026_06_12:
        def remap(ph):
            if ph == 1:
                return 10
            elif ph == 10:
                return 11
            elif ph == 11:
                return 1
            return ph

        iphsub = remap(iphsub)

        new_zones_age = []
        for ictherm, age_1, tp1, tp2, nph_layer, hcs, iphs, ixtb1, ixtb2 in zones_age:
            new_iphs = [remap(iph) for iph in iphs]
            new_zones_age.append((ictherm, age_1, tp1, tp2, nph_layer, hcs, new_iphs, ixtb1, ixtb2))
        zones_age = new_zones_age

        if len(phases) >= 11:
            new_phases = list(phases)
            new_phases[9] = phases[0]
            new_phases[10] = phases[9]
            new_phases[0] = phases[10]
            phases = new_phases

        new_inhoms = []
        for ix1_i, ix2_i, iy1_i, iy2_i, inphase_i, igeom_i, amp in inhoms:
            new_inhoms.append((ix1_i, ix2_i, iy1_i, iy2_i, remap(inphase_i), igeom_i, amp))
        inhoms = new_inhoms

    data = {
        'date': date,
        'nex': nex, 'nez': nez, 'x0': x0, 'z0': z0, 'rxbo': rxbo, 'rzbo': rzbo,
        'ircoord': ircoord, 'coordfile': coordfile, 'nzonx': nzonx, 'zones_x': zones_x,
        'nzony': nzony, 'zones_y': zones_y,
        'nystressbc': nystressbc, 'nydrsides': nydrsides, 'nofbc': nofbc, 'bcs': bcs,
        'nyhydro': nyhydro, 'pisos': pisos, 'iphsub': iphsub, 'drosub': drosub, 'damp_vis': damp_vis,
        'g': g, 'i_prestress': i_prestress, 'extra_pres': extra_pres, 'itherm': itherm,
        'istress_therm': istress_therm, 'ishearh': ishearh, 't_top': t_top, 't_bot': t_bot,
        'hs': hs, 'hr': hr, 'itemp_bc': itemp_bc, 'bot_bc': bot_bc, 'irtemp': irtemp,
        'tempfile': tempfile, 'nzone_age': nzone_age, 'zones_age': zones_age,
        'nphase': nphase, 'phases': phases, 'irphase': irphase, 'phasefile': phasefile,
        'inhom': inhom, 'inhoms': inhoms, 'ten_off': ten_off, 'tau_heal': tau_heal,
        'v_min': v_min, 'v_max': v_max, 'ivis_shape': ivis_shape, 'efoldc': efoldc,
        'itype_melting': itype_melting, 'nelem_serp': nelem_serp, 'nelem_dike': nelem_dike,
        'prod_magma': prod_magma, 'prod_magma2': prod_magma2, 'rho_magma': rho_magma, 'angle_mzone': angle_mzone,
        'fmagma_max': fmagma_max, 'ratio_mantle_mzone': ratio_mantle_mzone, 'ratio_crust_mzone': ratio_crust_mzone,
        'latent_heat_magma': latent_heat_magma, 'lambda_freeze': lambda_freeze,
        'lambda_freeze_tdep': lambda_freeze_tdep, 'weaken_ratio_plastic': weaken_ratio_plastic,
        'weaken_ratio_viscous': weaken_ratio_viscous, 'ny_rem': ny_rem, 'mode_rem': mode_rem,
        'ntest_rem': ntest_rem, 'angle_rem': angle_rem, 'dx_rem': dx_rem, 'topo_kappa': topo_kappa,
        'fac_kappa': fac_kappa, 'idt_scale': idt_scale, 'dt_scale': dt_scale, 'tolerance': tolerance,
        'i_rey': i_rey, 'xReyn': xReyn, 'ifreq_rmasses': if_rmasses, 'ifreq_imasses': if_imasses,
        'ifreq_visc': if_visc, 'ifreq_avgsr': if_avgsr, 'frac_elastic': frac_elastic,
        'frac_maxwell': frac_maxwell, 'movegrid': movegrid, 'ndim': ndim, 'damping': damping,
        'mix1': mix1, 'mix2': mix2, 'time_max': time_max, 'dtout_screen': dtout_screen,
        'dtout_file': dtout_file, 'io_flags': io_flags, 'lastout': lastout,
        'dtsave_file': dtsave_file, 'lastsave': lastsave
    }
    return data, comments_map

def write_new_inp(data, comments_map, final_comments, filename):
    with open(filename, 'w') as f:
        date = data['date']
        
        def remap_comment_text(comment):
            if not comment or data['nphase'] < 11:
                return comment
            # Swap (1) -> (10), (10) -> (11), (11) -> (1)
            comment = comment.replace('(1)', '__T1__')
            comment = comment.replace('(10)', '__T10__')
            comment = comment.replace('(11)', '__T11__')
            comment = comment.replace('__T1__', '(10)')
            comment = comment.replace('__T10__', '(11)')
            comment = comment.replace('__T11__', '(1)')
            return comment

        def write_line(key, content_str, default_preceding=None, default_trailing=None):
            if key in comments_map:
                preceding, trailing = comments_map[key]
                # Write preceding comments
                for is_blank, line in preceding:
                    if is_blank:
                        f.write('\n')
                    else:
                        line_text = line
                        if date < D_2026_06_12:
                            line_text = remap_comment_text(line_text)
                        f.write(line_text + '\n')
                # Write content and trailing
                trailing_text = trailing
                if date < D_2026_06_12:
                    trailing_text = remap_comment_text(trailing_text)
                if trailing_text and not trailing_text[0].isspace():
                    f.write(f"{content_str}  {trailing_text}\n")
                else:
                    f.write(f"{content_str}{trailing_text}\n")
            else:
                # Fallback to default
                if default_preceding:
                    for is_blank, line in default_preceding:
                        if is_blank:
                            f.write('\n')
                        else:
                            f.write(line + '\n')
                if default_trailing and not default_trailing[0].isspace():
                    f.write(f"{content_str}  {default_trailing}\n")
                else:
                    f.write(f"{content_str}{default_trailing or ''}\n")

        # 1. Mesh
        write_line('nex_nez', f"{data['nex']},{data['nez']}", 
                   default_preceding=[(False, ";=================================================================="),
                                     (False, ";             M e s h    P a r a m e t e r s "),
                                     (False, ";==================================================================")],
                   default_trailing="            number of _elements_ in X and Z directions")
        write_line('x0_z0', f"{data['x0']},{data['z0']}",
                   default_trailing="           x0,z0 begin.coord")
        write_line('rxbo_rzbo', f"{data['rxbo']},{data['rzbo']}",
                   default_trailing="    rxbo,rzbo (size of the region)")
        coord_val = f"{data['ircoord']}"
        if data['coordfile']:
            coord_val += f", {data['coordfile']}"
        write_line('ircoord', coord_val,
                   default_preceding=[(True, "")],
                   default_trailing="    ircoord, coordfile: read init. coordinates from file")
        
        write_line('nzonx', f"{data['nzonx']}",
                   default_preceding=[(True, ""), (False, "; X direction")],
                   default_trailing="     Number zones X-direction (0 - regular grid)")
        for idx in range(data['nzonx']):
            nelz, sizez = data['zones_x'][idx]
            write_line(f'zones_x_{idx}', f" {nelz}  {sizez}")
            
        write_line('nzony', f"{data['nzony']}",
                   default_preceding=[(True, ""), (False, "; Z direction")],
                   default_trailing="     Number zones Z-direction (0 - regular grid)")
        for idx in range(data['nzony']):
            nelz, sizez = data['zones_y'][idx]
            write_line(f'zones_y_{idx}', f" {nelz}  {sizez}")

        # 2. Mechanical boundary conditions
        write_line('nystressbc_nydrsides', f"{data['nystressbc']}  {data['nydrsides']}",
                   default_preceding=[(True, ""),
                                     (False, ";==================================================================="),
                                     (False, ";        C o n d i t i o n s:  M e c h a n i c a l"),
                                     (False, ";===================================================================")],
                   default_trailing="    nystressbc, nydrsides")
        write_line('nofbc', f"{data['nofbc']}",
                   default_trailing="    nofbc: Number of boundary conditions")
        for idx in range(data['nofbc']):
            bc = data['bcs'][idx]
            line_str = "  ".join(f"{v:g}" if isinstance(v, float) else str(v) for v in bc)
            write_line(f'bc_{idx}', line_str,
                       default_preceding=[(False, ";nofside  nbc1 nbc2  nbc   a       b    c     d     e     f      g     h     i ")] if idx == 0 else None)
                       
        write_line('nyhydro_pisos_iphsub_drosub_damp_vis',
                   f"{data['nyhydro']}  {data['pisos']}  {data['iphsub']}  {data['drosub']}  {data['damp_vis']}",
                   default_preceding=[(True, ""), (False, "; Hydrostatic pressure applied at bottom")])
                   
        write_line('g', f"{data['g']}",
                   default_preceding=[(True, ""), (False, "; Gravity")])

        # 3. Thermal
        write_line('i_prestress', f"{data['i_prestress']}",
                   default_preceding=[(True, ""),
                                     (False, ";============================================================="),
                                     (False, ";            C o n d i t i o n s : T h e r m a l"),
                                     (False, ";=============================================================")],
                   default_trailing="       -iprestress: allow topo build up by isostasy")
        write_line('extra_pres', f"{data['extra_pres']}",
                   default_trailing="    -extra pressure (Pa)")
        write_line('itherm', f"{data['itherm']}",
                   default_trailing="        -itherm (1-mech+therm, 2-no mech)")
        write_line('istress_therm', f"{data['istress_therm']}",
                   default_trailing="        -istress_therm: Add thermal stresses")
        write_line('ishearh', f"{data['ishearh']}",
                   default_trailing="        -ishearh: Add shear heating")
        write_line('t_top', f"{data['t_top']}",
                   default_trailing="       -t_top (Surface temperature)")
        write_line('t_bot', f"{data['t_bot']}",
                   default_trailing="     -t_bot (Bottom temperature)")
                   
        write_line('hs_hr', f"{data['hs']}, {data['hr']}",
                   default_preceding=[(True, ""), (False, "; Radiogenic heating")],
                   default_trailing="      -hs (W/kg), hr (km)")
                   
        write_line('itemp_bc_bot_bc', f"{data['itemp_bc']} {data['bot_bc']}",
                   default_preceding=[(True, ""), (False, "; Bottom Boundary condition flag and value")])
                   
        write_line('irtemp', f"{data['irtemp']}",
                   default_preceding=[(True, ""), (False, "; Predefined distributions")],
                   default_trailing="              irtemp (0,1)")
        write_line('tempfile', f"{data['tempfile']}")
        
        write_line('nzone_age', f"{data['nzone_age']}",
                   default_trailing="              - nzone_age")
        
        for idx in range(data['nzone_age']):
            ictherm, age_1, tp1, tp2, nph_layer, hcs, iphs, ixtb1, ixtb2 = data['zones_age'][idx]
            write_line(f'zone_age_{idx}_line1', f" {ictherm}, {age_1}, {tp1}, {tp2}, {ixtb1}, {ixtb2}")
            hc_line = " ".join(map(str, hcs[:nph_layer-1]))
            write_line(f'zone_age_{idx}_line2', f"    {nph_layer}, {hc_line}")
            iph_line = " ".join(map(str, iphs[:nph_layer]))
            write_line(f'zone_age_{idx}_line3', f"    {iph_line}")
            
        # 4. Rheology
        write_line('nphase', f"{data['nphase']}",
                   default_preceding=[(True, ""),
                                     (False, ";==================================================================="),
                                     (False, ";                     R h e o l o g y"),
                                     (False, ";===================================================================")],
                   default_trailing="  Number of Different Rheologies")
        
        for idx in range(data['nphase']):
            p = data['phases'][idx]
            line_str = ", ".join(f"{v:g}" if isinstance(v, float) else str(v) for v in p)
            
            comment_key = f'phase_{idx}'
            if date < D_2026_06_12 and data['nphase'] >= 11:
                if idx == 0:
                    comment_key = 'phase_10'
                elif idx == 9:
                    comment_key = 'phase_0'
                elif idx == 10:
                    comment_key = 'phase_9'
            
            write_line(comment_key, line_str,
                       default_preceding=[(False, ";irheol,_,den, alfa,  beta,    n,       A,       E,    V,    Lame:rl, Lame:rm,pls1,pls2,fric1,fric2, coh1, coh2,dilat1,dilat2,cond,    cp,     Ts,     Tl,     Tk, erodibility")] if idx == 0 else None)

        # 5. Phase distribution & Inhomogeneities
        write_line('irphase', f"{data['irphase']}",
                   default_preceding=[(True, ""), (False, "; INITIAL PHASE DISTRIBUTION")],
                   default_trailing="              ; irphase")
        write_line('phasefile', f"{data['phasefile']}")
        
        write_line('inhom', f"{data['inhom']}",
                   default_preceding=[(True, ""), (False, "; Initial heterogeneities")],
                   default_trailing="  - inhom(number of heterogeneities)")
        for idx in range(data['inhom']):
            inh = data['inhoms'][idx]
            line_str = "  ".join(f"{v:g}" if isinstance(v, float) else str(v) for v in inh)
            write_line(f'inhom_{idx}', f"  {line_str}")

        # 6. Mechanical limits & Healing
        write_line('ten_off', f"{data['ten_off']}",
                   default_preceding=[(True, ""), (False, "; Tension cut off")])
        write_line('tau_heal', f"{data['tau_heal']}",
                   default_preceding=[(True, ""), (False, "; linear healing parameter")])
        write_line('v_lims', f"{data['v_min']}, {data['v_max']}, {data['ivis_shape']}, {data['efoldc']}",
                   default_preceding=[(True, ""), (False, "; VISCOSITY LIMIT")])
                   
        write_line('magma_line1', f"{data['itype_melting']}, {data['nelem_serp']}, {data['nelem_dike']}, {data['prod_magma']}, {data['prod_magma2']}, {data['rho_magma']}",
                   default_preceding=[(True, ""), (False, "; Magma:"), (False, "; itype_melting, nelem_serp, nelem_dike, prod_magma, prod_magma2, rho_magma")])
        write_line('magma_line2', f"{data['angle_mzone']}, {data['fmagma_max']}, {data['ratio_mantle_mzone']}, {data['ratio_crust_mzone']}",
                   default_preceding=[(False, "; angle_mzone, fmagma_max, ratio_mantle_mzone, ratio_crust_mzone")])
        write_line('magma_line3', f"{data['latent_heat_magma']}, {data['lambda_freeze']}, {data['lambda_freeze_tdep']}",
                   default_preceding=[(False, "; latent_heat_magma, lambda_freeze, lambda_freeze_tdep")])
        write_line('magma_line4', f"{data['weaken_ratio_plastic']}, {data['weaken_ratio_viscous']}",
                   default_preceding=[(False, "; weaken_ratio_plastic, weaken_ratio_viscous")])

        # 7. Remeshing
        write_line('rem_line', f"{data['ny_rem']}  {data['mode_rem']}  {data['ntest_rem']}  {data['angle_rem']}",
                   default_preceding=[(True, ""),
                                     (False, ";================================================================="),
                                     (False, ";                       R e m e s h i n g"),
                                     (False, ";=================================================================")])
        write_line('dx_rem', f"{data['dx_rem']}",
                   default_preceding=[(False, "; dx_rem")])
        write_line('topo_kappa_fac_kappa', f"{data['topo_kappa']}  {data['fac_kappa']}",
                   default_preceding=[(False, "; topo diffusivity and amplification factor")])

        # 8. Process control
        write_line('idt_scale', f"{data['idt_scale']}",
                   default_preceding=[(True, ""),
                                     (False, ";================================================================="),
                                     (False, ";                   P r o c e s s   c o n t r o l"),
                                     (False, ";=================================================================")],
                   default_trailing="         - idt_scale")
        write_line('dt_scale_tolerance', f"{data['dt_scale']},{data['tolerance']}",
                   default_trailing="  - dt_scale, tolerance")
        write_line('i_rey_xReyn', f"{data['i_rey']},{data['xReyn']}",
                   default_trailing="   - i_rey, Reyn")
        write_line('if_rmasses', f"{data['ifreq_rmasses']}",
                   default_trailing="       - frequency of evaluation of real masses")
        write_line('if_imasses', f"{data['ifreq_imasses']}",
                   default_trailing="       - frequency of evaluation of inertial masses")
        write_line('if_visc', f"{data['ifreq_visc']}",
                   default_trailing="       - frequency of evaluation of Non-Newtonian visc")
        write_line('if_avgsr', f"{data['ifreq_avgsr']}",
                   default_trailing="       - frequency of averaging strain rate")
                   
        write_line('frac_elastic_frac_maxwell', f"{data['frac_elastic']},{data['frac_maxwell']}",
                   default_preceding=[(True, "")])
        write_line('movegrid_ndim', f"{data['movegrid']},{data['ndim']}")
        write_line('damping_mix1_mix2', f"{data['damping']},{data['mix1']},{data['mix2']}")

        # 9. Output
        write_line('time_max', f"{data['time_max']}",
                   default_preceding=[(True, ""),
                                     (False, ";======================================================================\n"),
                                     (False, ";                             O U T P U T\n"),
                                     (False, ";======================================================================\n")],
                   default_trailing="  ; max time in kyrs")
        write_line('dtout_screen', f"{data['dtout_screen']}",
                   default_trailing="  ; dtout_screen in kyrs")
        write_line('dtout_file', f"{data['dtout_file']}",
                   default_trailing="  ; dtout_file in kyrs")
                   
        io_str = "   ".join(map(str, data['io_flags']))
        write_line('io_flags', f"  {io_str}",
                   default_preceding=[(True, ""), (False, "; Variables to print")])
                   
        write_line('lastout', f"{data['lastout']}",
                   default_preceding=[(True, ""), (False, "; output for last step only (1) or each nout step (0)")])
        write_line('dtsave_file', f"{data['dtsave_file']}",
                   default_preceding=[(True, ""), (False, "; Time interval for process saving")])
        write_line('lastsave', f"{data['lastsave']}",
                   default_preceding=[(True, ""), (False, "; saving last step only (1) or each nsave step (0)")])

        for is_blank, line in final_comments:
            if is_blank:
                f.write('\n')
            else:
                line_text = line
                if date < D_2026_06_12:
                    line_text = remap_comment_text(line_text)
                f.write(line_text + '\n')

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 update_inp.py <input_filename> <YYYY-MM-DD_code_date>")
        sys.exit(1)
        
    filename = sys.argv[1]
    date_str = sys.argv[2]
    
    try:
        date = parse_date(date_str)
    except ValueError as e:
        print(e)
        sys.exit(1)
        
    print(f"Tokenizing and parsing '{filename}' targeting code date {date.strftime('%Y-%m-%d')}...")
    clean_lines, final_comments = get_clean_lines(filename)
    
    try:
        data, comments_map = parse_inp(clean_lines, date)
    except Exception as e:
        print(f"Error parsing input file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
        
    # Write to a new file named with '_upgraded.inp' suffix
    out_filename = filename
    if filename.endswith('.inp'):
        out_filename = filename[:-4] + '_upgraded.inp'
    else:
        out_filename = filename + '_upgraded.inp'
        
    print(f"Writing upgraded input file to '{out_filename}'...")
    write_new_inp(data, comments_map, final_comments, out_filename)
    print("Upgrade completed successfully!")

if __name__ == '__main__':
    main()
