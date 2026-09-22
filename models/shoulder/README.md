# solutionspace

model: das3_PhD.osim is renamed from das3_PhD_v1.osim from Study1 folder to keep filename consistent
- Update (2024-05-21):
    - Changed the wrapping object for the rotator cuff muscles.  Added a sphere (based on MSc work) called "rotcuff" for each of the RC muscles.
    - Added the "mid-point articular surface glenoid" as Glenoid marker -- to be used for orienting the stability ratios to glenoid plane.
        - Was based on l1091.dsp file (coordinates correspond to global in osim --> then converted to scapula_r using osim gui)
