# References

- `csr_t14_5_02.csv`: copied by hand from the CDISC Pilot CSR (`cdiscpilot01.pdf`, PDF page 79,
  Table 14-5.02, "Incidence of Treatment Emergent Serious Adverse Events by Treatment Group").
  Denominators: Placebo N=86, Xanomeline Low N=84, Xanomeline High N=84. Blank cells in the CSR
  are blank here. `trtan` uses the ADaM codes (0, 54, 81).
- The `adae` reference is the released ADaM `adae.xpt`, downloaded by `bench/fetch.py`.
- `csr_t14_2_01.csv`: CSR Table 14-2.01 (PDF pages 46–48), parsed from the CSR's text layout and
  reconciled by count; same file as `forward2028/tlf/refs/t14_2_01.csv`.
