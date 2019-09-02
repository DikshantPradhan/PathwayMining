\* Problem: FBA_test_model *\

Maximize
 obj: 0 v4

Subject To
 S: - v3 - v2 + v1 = 0
 P: + v4 - 2 v3 - v2 = 0
 S_P: - v5 + v2 = 0
 S_P_P: - v6 + v3 = 0

Bounds
 0 <= v1 <= 1000
 0 <= v2 <= 1000
 0 <= v3 <= 1000
 0 <= v4 <= 1000
 0 <= v5 <= 1000
 0 <= v6 <= 1000

End
