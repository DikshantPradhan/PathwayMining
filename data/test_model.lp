\* Problem: FBA_test_model *\

Maximize
 obj: 0 v1

Subject To
 S: + v1 - v2 - v3 = 0
 P: + v4 - v2 - 2 v3 = 0
 S_P: + v2 - v5 = 0
 S_P_P: + v3 - v6 = 0

Bounds
 0 <= v1 <= 1000
 -1000 <= v2 <= 1000
 -1000 <= v3 <= 1000
 0 <= v4 <= 1000
 -1000 <= v5 <= 1000
 -1000 <= v6 <= 1000


End
