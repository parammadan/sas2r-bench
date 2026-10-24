# References

- `csr_t14_5_02.csv`: copied by hand from the CDISC Pilot CSR (`cdiscpilot01.pdf`, PDF page 79,
  Table 14-5.02, "Incidence of Treatment Emergent Serious Adverse Events by Treatment Group").
  Denominators: Placebo N=86, Xanomeline Low N=84, Xanomeline High N=84. Blank cells in the CSR
  are blank here. `trtan` uses the ADaM codes (0, 54, 81).
- The `adae` reference is the released ADaM `adae.xpt`, downloaded by `bench/fetch.py`.
- `csr_t14_2_01.csv`: CSR Table 14-2.01 (PDF pages 46–48), parsed from the CSR's text layout and
  reconciled by count; same file as `forward2028/tlf/refs/t14_2_01.csv`.
- `csr_t14_5_01.csv`: CSR Table 14-5.01 (PDF pages 63–78), parsed from the CSR's text layout.
  Wrapped names are joined; a row is a system organ class when it equals an ADAE AEBODSYS value
  or is one cut at the 40-character column. Parsing QA against the released ADAE: 685 of 690
  preferred-term counts equal it directly; the other 5 are names the CSR cut at 38 characters.
- `csr_t14_1_01.csv`: CSR Table 14-1.01 (PDF page 43), same file as `forward2028/tlf/refs/t14_1_01.csv`.
- `expected/<id>.json`: datasets each structure-checked program creates, from `bench/expected.py`.
