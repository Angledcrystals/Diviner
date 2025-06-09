#!/usr/bin/env python3
"""
G-S Alignment Data Converter
Converts G-S alignment data dump to JSON format for Divine Pixel Tool

Author: Angledcrystals
Date: 2025-06-09
Compatible with: G-S Stereo Viewer & Divine Pixel Tool v1.2
"""

import json
import re
import sys
from datetime import datetime

def parse_gs_alignment_dump(file_path):
    """Parse G-S alignment data dump and convert to JSON format."""
    
    print(f"🔄 Converting G-S alignment data: {file_path}")
    
    alignments = []
    current_condition = None
    header_parsed = False
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            
            # Skip empty lines and header info
            if not line or line.startswith('S-NUIT') or line.startswith('=') or line.startswith('Generated:'):
                continue
            
            # Parse condition headers
            if line.startswith('CONDITION:'):
                current_condition = line.replace('CONDITION:', '').strip()
                header_parsed = False
                print(f"   📊 Processing condition: {current_condition}")
                continue
            
            # Look for CSV header line
            if 'G_theta_deg,G_phi_deg,Hadit_theta_deg,Hadit_phi_deg,S_x,S_y' in line:
                header_parsed = True
                continue
            
            # Skip summary and statistical lines
            if any(keyword in line for keyword in [
                'Number of alignments:', 'STATISTICAL SUMMARY:', 'Mean:', 'Min:', 'Max:',
                'range:', 'BEST ALIGNMENT:', 'DETAILED ALIGNMENT DATA:', '----'
            ]):
                continue
            
            # Parse alignment data lines (CSV format)
            if header_parsed and current_condition and ',' in line:
                try:
                    # Split CSV line
                    parts = line.split(',')
                    
                    # Ensure we have enough parts
                    if len(parts) >= 6:
                        # Extract the required fields (first 6 columns)
                        g_theta = float(parts[0])
                        g_phi = float(parts[1])
                        hadit_theta = float(parts[2])
                        hadit_phi = float(parts[3])
                        s_x = float(parts[4])
                        s_y = float(parts[5])
                        
                        # Create alignment object in G-S format
                        alignment = {
                            'S_x': s_x,
                            'S_y': s_y,
                            'G_theta': g_theta,
                            'G_phi': g_phi,
                            'Hadit_theta': hadit_theta,
                            'Hadit_phi': hadit_phi,
                            'condition': current_condition,  # Optional metadata
                        }
                        
                        # Add optional fields if available
                        if len(parts) >= 7:
                            alignment['dist_boundary'] = float(parts[6])
                        if len(parts) >= 8:
                            alignment['dist_origin'] = float(parts[7])
                        if len(parts) >= 11:
                            alignment['G_refl_x'] = float(parts[8])
                            alignment['G_refl_y'] = float(parts[9])
                            alignment['G_refl_z'] = float(parts[10])
                        
                        alignments.append(alignment)
                        
                except (ValueError, IndexError) as e:
                    print(f"   ⚠️ Warning: Could not parse line {line_num}: {line[:50]}...")
                    print(f"      Error: {e}")
                    continue
    
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return None
    
    print(f"✅ Successfully converted {len(alignments)} alignments")
    return alignments

def save_json_alignment_data(alignments, output_path):
    """Save alignments in JSON format."""
    
    try:
        # Create metadata
        metadata = {
            'generated_by': 'G-S Alignment Data Converter',
            'author': 'Angledcrystals',
            'conversion_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC'),
            'total_alignments': len(alignments),
            'format_version': '1.0',
            'compatible_with': ['G-S Stereo Viewer', 'Divine Pixel Tool v1.2'],
            'required_fields': ['S_x', 'S_y', 'G_theta', 'G_phi'],
            'optional_fields': ['Hadit_theta', 'Hadit_phi', 'condition', 'dist_boundary', 'dist_origin']
        }
        
        # Create output structure
        output_data = {
            'metadata': metadata,
            'alignments': alignments
        }
        
        # Save with pretty formatting
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(output_data, f, indent=2, ensure_ascii=False)
        
        print(f"💾 Saved JSON alignment data: {output_path}")
        return True
        
    except Exception as e:
        print(f"❌ Error saving JSON: {e}")
        return False

def save_simple_json_alignment_data(alignments, output_path):
    """Save alignments in simple JSON array format (for direct tool use)."""
    
    try:
        # Save as simple array (exactly what Divine Pixel Tool expects)
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(alignments, f, indent=2, ensure_ascii=False)
        
        print(f"💾 Saved simple JSON alignment data: {output_path}")
        return True
        
    except Exception as e:
        print(f"❌ Error saving simple JSON: {e}")
        return False

def analyze_alignment_data(alignments):
    """Analyze the converted alignment data."""
    
    if not alignments:
        return
    
    print("\n📊 G-S ALIGNMENT DATA ANALYSIS")
    print("=" * 50)
    
    # Basic statistics
    print(f"Total alignments: {len(alignments)}")
    
    # Field analysis
    s_x_values = [a['S_x'] for a in alignments]
    s_y_values = [a['S_y'] for a in alignments]
    g_theta_values = [a['G_theta'] for a in alignments]
    g_phi_values = [a['G_phi'] for a in alignments]
    
    print(f"S_x range: [{min(s_x_values):.6f}, {max(s_x_values):.6f}]")
    print(f"S_y range: [{min(s_y_values):.6f}, {max(s_y_values):.6f}]")
    print(f"G_theta range: [{min(g_theta_values):.1f}°, {max(g_theta_values):.1f}°]")
    print(f"G_phi range: [{min(g_phi_values):.1f}°, {max(g_phi_values):.1f}°]")
    
    # Check for Hadit data
    hadit_count = sum(1 for a in alignments if 'Hadit_theta' in a and 'Hadit_phi' in a)
    print(f"Alignments with Hadit data: {hadit_count}/{len(alignments)}")
    
    if hadit_count > 0:
        hadit_theta_values = [a['Hadit_theta'] for a in alignments if 'Hadit_theta' in a]
        hadit_phi_values = [a['Hadit_phi'] for a in alignments if 'Hadit_phi' in a]
        print(f"Hadit_theta range: [{min(hadit_theta_values):.1f}°, {max(hadit_theta_values):.1f}°]")
        print(f"Hadit_phi range: [{min(hadit_phi_values):.1f}°, {max(hadit_phi_values):.1f}°]")
    
    # Condition analysis
    conditions = set(a.get('condition', 'Unknown') for a in alignments)
    print(f"Conditions found: {len(conditions)}")
    for condition in sorted(conditions):
        count = sum(1 for a in alignments if a.get('condition') == condition)
        print(f"  • {condition}: {count} alignments")
    
    # Sample alignment
    print(f"\nSample alignment:")
    sample = alignments[0]
    for key, value in sample.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.6f}")
        else:
            print(f"  {key}: {value}")

def main():
    """Main conversion function."""
    
    print("🔮 G-S Alignment Data Converter")
    print("Converting alignment data dump to JSON format for Divine Pixel Tool")
    print("Author: Angledcrystals")
    print("Date: 2025-06-09 09:00:42 UTC")
    print("=" * 70)
    
    # Get input file
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = input("📁 Enter path to alignment data dump file: ").strip()
    
    if not input_file:
        input_file = "alignmentdatadump.txt"  # Default filename
    
    # Parse alignment data
    alignments = parse_gs_alignment_dump(input_file)
    
    if not alignments:
        print("❌ No alignment data found or conversion failed")
        return
    
    # Analyze data
    analyze_alignment_data(alignments)
    
    # Generate output filenames
    base_name = input_file.rsplit('.', 1)[0]
    json_file = f"{base_name}_gs_alignments.json"
    simple_json_file = f"{base_name}_simple.json"
    
    # Save both formats
    print(f"\n💾 SAVING CONVERTED DATA")
    print("-" * 30)
    
    # Save full format with metadata
    if save_json_alignment_data(alignments, json_file):
        print(f"✅ Full format: {json_file}")
    
    # Save simple format for direct tool use
    if save_simple_json_alignment_data(alignments, simple_json_file):
        print(f"✅ Simple format: {simple_json_file}")
    
    print(f"\n🔮 READY FOR DIVINE PIXEL TOOL")
    print("=" * 40)
    print(f"Use this file in Divine Pixel Tool: {simple_json_file}")
    print("The tool now supports your exact G-S alignment data format!")
    print("\n📋 Next steps:")
    print("1. Open Divine Pixel Tool")
    print("2. Load your image to divine")
    print(f"3. Load G-S alignment data: {simple_json_file}")
    print("4. Enable 'Use G-S Alignment Data'")
    print("5. Generate G-S coordinate maps")
    print("6. Choose 'gs_enhanced_divine' method")
    print("7. Divine your pixels with G-S enhancement!")

if __name__ == "__main__":
    main()
