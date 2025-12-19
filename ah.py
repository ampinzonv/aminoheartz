import numpy as np
from scipy.io import wavfile
import argparse
import os
import matplotlib.pyplot as plt
from midiutil import MIDIFile 

def read_fasta(file_path):
    """Reads a FASTA file and returns the clean amino acid sequence."""
    sequence = ""
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith(">"):
                    sequence += line.strip()
        return sequence.upper()
    except FileNotFoundError:
        print(f"Error: The file {file_path} does not exist.")
        return None

def get_window_counts(window_seq):
    """Calculates the relative intensity of each bond in a specific window."""
    keys = ["Amide_I (C=O)", "Amide_II (N-H)", "C-H Aliphatic", "O-H/N-H (Polar)", "Disulfide Bridge", "Ionic (Charged)"]
    counts = {k: 0.0 for k in keys}
    w_len = max(len(window_seq), 1)
    
    for aa in window_seq:
        counts["Amide_I (C=O)"] += 1
        counts["Amide_II (N-H)"] += 1
        
        if aa in "AVILMF": counts["C-H Aliphatic"] += 2
        if aa in "STNQYC": counts["O-H/N-H (Polar)"] += 1
        if aa in "C": 
            counts["Disulfide Bridge"] += 0.5
            counts["O-H/N-H (Polar)"] += 1
        if aa in "DEKRH":
            counts["Ionic (Charged)"] += 2

    return {k: v / w_len for k, v in counts.items()}

def print_sequence_stats(control_points):
    """Prints a formatted table of bond statistics and intensities."""
    # Metadatos para la tabla (Frecuencia y Nota Musical)
    # Nota: Se añade Ionic (Charged) que faltaba en el ejemplo manual
    bond_meta = {
        "Amide_I (C=O)":       (330.0, "E 4"),
        "Amide_II (N-H)":      (311.1, "D# 4"),
        "C-H Aliphatic":       (587.3, "D 5"),
        "O-H/N-H (Polar)":     (659.2, "E 5"),
        "Disulfide Bridge":    (123.4, "B 2"),
        "Ionic (Charged)":     (1400.0, "F 6")
    }

    # Calcular promedios globales
    num_points = len(control_points)
    if num_points == 0: return

    totals = {k: 0.0 for k in bond_meta.keys()}
    for point in control_points:
        for k, v in point.items():
            totals[k] += v

    print("\n" + "="*65)
    print(f"{'BOND / VIBRATION':<20} | {'FREQ (Hz)':<9} | {'NOTE':<6} | {'INTENSITY'}")
    print("-" * 65)

    for bond, (freq, note) in bond_meta.items():
        # Promedio de intensidad en toda la proteína
        avg_intensity = totals[bond] / num_points
        
        # Crear barra visual (escala x20 para que se vea bien en consola)
        bar_len = int(avg_intensity * 20)
        bar = "█" * bar_len
        
        print(f"{bond:<20} | {freq:<9.1f} | {note:<6} | {avg_intensity:.3f} {bar}")
    
    print("="*65 + "\n")

def plot_protein_profile(control_points, output_name):
    """Generates a PNG graph of the chemical intensities."""
    plt.figure(figsize=(12, 6))
    bonds = control_points[0].keys()
    x = np.arange(len(control_points))
    
    for bond in bonds:
        y = [point[bond] for point in control_points]
        plt.plot(x, y, label=bond, linewidth=1.5, alpha=0.8)
    
    plt.title(f"AminoHeartz: Bio-Physical Profile of {output_name}")
    plt.xlabel("Window Step")
    plt.ylabel("Intensity")
    plt.legend(loc='upper right', fontsize='small')
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    graph_filename = output_name.replace(".wav", ".png").replace(".mid", ".png")
    plt.savefig(graph_filename, dpi=150)
    print(f"📊 Visualization saved as: {graph_filename}")

def generate_midi_file(control_points, total_duration, filename):
    """Generates a MIDI file with sustained notes and CC11 Expression curves."""
    midi = MIDIFile(6)
    
    midi_map = {
        "Amide_I (C=O)":       (64, 0), 
        "Amide_II (N-H)":      (63, 1),
        "C-H Aliphatic":       (74, 2),
        "O-H/N-H (Polar)":     (76, 3),
        "Disulfide Bridge":    (47, 4),
        "Ionic (Charged)":     (89, 5)
    }

    tempo = 60 

    for bond, (note_num, track) in midi_map.items():
        midi.addTrackName(track, 0, bond)
        midi.addTempo(track, 0, tempo)
        
        midi.addNote(track, 0, note_num, 0, total_duration, 90)

        times = np.linspace(0, total_duration, len(control_points))
        for i, t in enumerate(times):
            intensity = control_points[i][bond]
            cc_val = int(min(intensity * 60, 127))
            midi.addControllerEvent(track, 0, float(t), 11, cc_val)

    midi_filename = filename.replace(".wav", ".mid")
    with open(midi_filename, "wb") as output_file:
        midi.writeFile(output_file)
    print(f"🎹 MIDI file generated: {midi_filename}")

def generate_evolving_sound(sequence, filename, window_size=10, speed_aa_per_sec=4.0, make_midi=False, show_stats=False):
    sample_rate = 44100
    total_duration = len(sequence) / speed_aa_per_sec
    total_samples = int(sample_rate * total_duration)
    t = np.linspace(0, total_duration, total_samples)
    
    print(f"\n--- 🧬 Starting AminoHeartz Scanner ---")
    print(f"Sequence Length: {len(sequence)} aa | Duration: {total_duration:.2f}s")
    
    bond_data = {
        "Amide_I (C=O)":       330.0,
        "Amide_II (N-H)":      311.1,
        "C-H Aliphatic":       587.3,
        "O-H/N-H (Polar)":     659.2,
        "Disulfide Bridge":    123.4,
        "Ionic (Charged)":     1400.0    
    }

    # 1. Calculate Control Points
    control_points = []
    for i in range(len(sequence)):
        start = i
        end = min(i + window_size, len(sequence))
        intensities = get_window_counts(sequence[start:end])
        control_points.append(intensities)

    # --- NUEVA FUNCIONALIDAD: ESTADÍSTICAS ---
    if show_stats:
        print_sequence_stats(control_points)
    # -----------------------------------------

    # 2. Visualization
    plot_protein_profile(control_points, filename)

    # 3. MIDI Generation
    if make_midi:
        generate_midi_file(control_points, total_duration, filename)

    # 4. Audio Synthesis
    final_audio = np.zeros(total_samples)
    x_control = np.linspace(0, total_duration, len(control_points))
    
    print("Synthesizing audio textures...")
    for bond, freq in bond_data.items():
        y_control = [point[bond] for point in control_points]
        envelope = np.interp(t, x_control, y_control)
        oscillator = np.sin(2 * np.pi * freq * t) * envelope
        final_audio += oscillator

    if np.max(np.abs(final_audio)) > 0:
        final_audio = final_audio / np.max(np.abs(final_audio))
        
    audio_int16 = (final_audio * 32767).astype(np.int16)
    wavfile.write(filename, sample_rate, audio_int16)
    print(f"✅ Audio file generated: {filename}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AminoHeartz: Protein Sonification")
    parser.add_argument("input", help="Path to input .fasta file")
    parser.add_argument("-o", "--output", help="Output .wav filename", default="protein_sound.wav")
    parser.add_argument("-w", "--window", type=int, help="Window size", default=10)
    parser.add_argument("-s", "--speed", type=float, help="Speed (aa/s)", default=5.0)
    parser.add_argument("-m", "--midi", action="store_true", help="Generate a .mid file for DAWs")
    parser.add_argument("--stats", action="store_true", help="Display biophysical statistics table")

    args = parser.parse_args()
    seq = read_fasta(args.input)
    if seq:
        generate_evolving_sound(seq, args.output, args.window, args.speed, args.midi, args.stats)