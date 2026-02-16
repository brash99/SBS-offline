# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SBS-offline is the reconstruction and analysis code for SuperBigBite (SBS) experiments at Jefferson Lab. It builds a single shared library (`libsbs.so`) containing detector classes, hardware decoders, and simulation support that plug into the Hall A Analyzer (Podd) framework.

## Build Commands

**Prerequisites:** ROOT 6+, Podd 1.7+ (pointed to by `$ANALYZER`), `$ROOTSYS` set.

**Full clean build** (from the SBS-offline source directory):
```bash
cd ../build && rm -fr * && cmake -DCMAKE_INSTALL_PREFIX=../install ../SBS-offline && make install -j8 && cd ../SBS-offline
```

**Incremental rebuild** (after code changes):
```bash
cd ../build && make install -j8 && cd ../SBS-offline
```

The build scripts `cmake_from_total_scratch` and `cmake_update_library` automate these with environment setup.

**CMake options:** `-DVERBOSE=ON` (default), `-DTESTCODE=ON` (default), `-DMCDATA=ON` (default), `-DGEN=true` (use GEN-replay instead of SBS-replay).

There is no formal test suite.

## Architecture

### Framework Integration

All detector and spectrometer classes inherit from Podd (Hall A Analyzer) base classes (`THaSpectrometer`, `THaDetector`, etc.). Podd handles the event loop, database loading, and output tree generation. SBS-offline registers its classes via ROOT dictionary generation (`sbs_LinkDef.h` → `sbsDict.cxx`).

### Key Class Hierarchies

**Spectrometers** — Top-level reconstruction orchestrators:
- `SBSBigBite` — Main BigBite spectrometer; coordinates tracking + PID
- `SBSEArm` / `SBSGEPEArm` — Electron arm reconstruction

**GEM Tracking** (the largest subsystem by code volume):
- `SBSGEMTrackerBase` — Core tracking algorithms (hit clustering, track finding/fitting)
- `SBSGEMSpectrometerTracker` / `SBSGEMPolarimeterTracker` — Concrete tracker implementations
- `SBSGEMModule` — Single GEM module readout and hit reconstruction (~44K lines)

**Generic Detector Framework:**
- `SBSGenericDetector` — Flexible base for detectors with ADC/TDC data; handles database-driven configuration of channel mapping, pedestals, gains
- `SBSElement` — Individual detector channel with associated `SBSData` (ADC/TDC/Waveform)
- Specialized subclasses: `SBSCalorimeter`, `SBSTimingHodoscope`, `SBSCherenkovDetector`, `SBSScintPlane`, `SBSCDet`

**Hardware Decoders** (in `Decoder::` namespace):
- `MPDModule` / `MPDModuleVMEv4` — Multi-Purpose Digitizer (GEM readout)
- `SBSDecodeF1TDCModule` — F1 TDC
- `VETROCModule`, `vfTDCModule`, `VTPModule` — Other digitizer types

**Simulation:**
- `SBSSimDecoder` — Converts Geant4 (g4sbs) digitized output into analyzer-compatible events
- `gmn_tree_digitized`, `genrp_tree_digitized`, `gep_tree_digitized` — Experiment-specific sim tree readers

### Data Flow

1. Raw data or simulation files → Hardware decoder modules → Detector `Decode()` methods
2. Detectors process hits via `CoarseProcess()` then `FineProcess()`
3. Spectrometer classes combine detector outputs for track reconstruction and PID
4. All detectors read calibration/geometry from TDatime-aware database files via `ReadDatabase()`

### Adding a New Detector

1. Create `SBSNewDet.h/.cxx` inheriting from `SBSGenericDetector` (or appropriate base)
2. Add the `.cxx` file to the `sources` list in `CMakeLists.txt`
3. Add `#pragma link C++ class SBSNewDet+;` to `sbs_LinkDef.h`
4. Implement `ReadDatabase()`, `Decode()`, `CoarseProcess()`, `FineProcess()` as needed

### Build Output Layout

After `make install`, the install directory contains:
- `lib/` — `libsbs.so` shared library
- `include/` — All headers + ROOT dictionary files
- `bin/` — `sbsenv.sh`/`sbsenv.csh` environment scripts
- `etc/` — `rootlogon.C` (auto-loads libsbs)
- `run_replay_here/` — `.rootrc` with macro paths configured

### Important Files

- `SBSData.h` — Core data structures (`ADCData`, `TDCData`, `WaveformData`, `PulseADCData`)
- `Helper.h` — Template utilities for container cleanup and vector initialization
- `DebugDef.h` — Debug level macros (0=none through 4=verbose)
- `sbs_LinkDef.h` — ROOT dictionary class registration (must be updated when adding classes)