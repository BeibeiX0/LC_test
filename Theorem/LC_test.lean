import Mathlib
import Theorem.valid.Defs_Combinatorial_Identities

open BigOperators Nat Finset Real

theorem test_001_1 {n : ℕ} : ∑ k ∈ Finset.range (n + 1), ((-1) ^ k * k * stirling_first n k : ℝ) = (-1 : ℝ)^n * (2^n - 1) := sorry

theorem test_002_1 {n : ℕ} : ∑ k ∈ Finset.range (n + 1),
  ((Nat.choose (2 * n) (n + k)) ^ 2 / (Nat.choose (2 * n + 2 * k) (n + k)) ^ 2 : ℝ) *
    ((4 * k + 1) * 2 ^ (4 * k) / (2 * n + 2 * k + 1) ^ 2 : ℝ) = 1 / (4 * n + 1) := sorry

theorem test_003_1 {n r : ℕ} : ∑ k ∈ Finset.range (n + 1), Nat.choose n k * Nat.choose k r * 2 ^ k = 2 ^ r * 3 ^ (n - r) * Nat.choose n r := sorry

theorem test_004_1 {x : ℝ} {n : ℕ} : ∑ k ∈ Finset.range (n + 1), ((-1) ^ k * Nat.choose n k * Nat.choose (2 * k) (n + k) * x ^ k : ℝ) = 2 ^ n * Nat.choose (2 * n) n := sorry

theorem test_005_1 {n : ℕ} : ∑ k ∈ Finset.range (n + 1),((-1) ^ k * Nat.choose n k  * k * 2 ^ k : ℝ) = 2 * n * (-1 : ℝ) ^ n := sorry

theorem test_006_1 {n : ℕ} :
    Real.log 2 = ∑' k : ℕ, (-1 : ℝ)^(k - 1) / k.succ := sorry

theorem test_007_1  {x : ℝ} (hx : |x| < 1) :
    Real.log (1 + x) = ∑' k : ℕ, (-1 : ℝ)^k * k * x^(k + 1) / (k + 1) := sorry

theorem test_008_1 (x y : ℕ):
    Complex.betaIntegral x y = 1 / (∑ k ∈ Finset.range (x + y), stirling_second (x + y - 1) k) := sorry

theorem test_009_1 {n : ℕ} (hn :n ≥ 1):stirling_second n 2=2^(n-1)-1 := sorry

theorem test_010_1 {n : ℕ} : ∑ k ∈ Finset.range (n + 1), k * (k !) = (n+1)! - 1 := sorry

theorem test_011_1 {n : ℕ} : ∑ k ∈ Finset.range (n + 1),
  (-1 : ℝ) ^ (n - k) * 2 ^ (2 * k) * Nat.choose (n + k + 1) (2 * k + 1)= n + 1 := sorry

theorem test_012_1 {n : ℕ} :
    catalan n = choose (2 * n) n - choose (2 * n) (n - 1) := sorry

theorem test_013_1 (n : ℕ) :
  ∑ k ∈ Finset.range (n + 1), ((k + 1) ^ 2 * Nat.choose n k : ℝ) = (2 : ℝ) ^ (n - 2) * (n ^ 2 + 5 * n + 4) := sorry

theorem test_014_1 (n:ℕ):∑ k ∈ Finset.range (n + 1), 1 / ((k + 1) * (k + 3)) * Nat.choose n k=(2^(n+3)-n-4)/(2*(n+2)*(n+3)) := sorry

theorem test_015_1 (n : ℕ) (h : n > 0) :
    ∑ k ∈ Finset.range (n + 1), Nat.choose (2 * n) (2 * k + 1) = 2 ^ (2 * n - 1) := sorry

theorem test_016_1 (n : ℕ) (x : ℝ) :
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * (n.choose k) * 2 ^ k * (cos (x / 2)) ^ k * cos ((k / 2) * x) =
      (-1 : ℝ) ^ n * cos (n * x) := sorry

theorem test_017_1 (n : ℕ) (hn : 4 ≤ n) :
    stirling_second n (n - 2) = n.choose 3 + 3 * n.choose 4 := sorry

theorem test_018_1 (n : ℕ) (hn : 1 ≤ n):
  ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * (k + 1) ^ n * Nat.choose (n + 1) (k + 1) = 0 := sorry

theorem test_019_1 (n m :ℕ):
  ∑ k ∈ Finset.range (n + 1), n * choose n k * choose m k = n * choose (n + m) n := sorry

theorem test_020_1 (n m r : ℕ) :
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * choose n k * choose (m - k) r = choose (m - n) (r - n) := sorry

theorem test_021_1 (n r : ℕ) (hn : 0 < n) (hrn : r ≤ n - 1) :
  Nat.choose n r = (n : ℝ) / (n - r) * Nat.choose (n - 1) r := sorry

theorem test_022_1 (n : ℕ):
  ∑ k ∈ Finset.range (n + 1), 1 / (((k : ℝ) + 1) * (k + 2)) * Nat.choose n k
  = (2 ^ (n + 2) - n - 3) / (((n : ℝ) + 1) * (n + 2)) := sorry

theorem test_023_1 (n : ℕ) :
  ∑ k ∈ Finset.range (n + 1),((k : ℝ) + 2) / (k + 1) * Nat.choose n k
  = (((n : ℝ) + 3) * 2 ^ n - 1) / (n + 1) := sorry

theorem test_024_1 (n : ℕ) :
  ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * 1 / (n+k+1) * Nat.choose n k = (n !) ^ 2 / (2 * n + 1) ! := sorry

theorem test_025_1(n:ℕ) : ∑ k ∈ Finset.Ico 1 (n+1),k^2*Nat.choose n k=n*(n+1)*2^(n-2) := sorry

theorem test_026_1(n m:ℕ)(hnm:m≤n):
  Nat.choose (n+r-1) r=∑ k ∈ Finset.range (m+1),Nat.choose (m+k+1) k *Nat.choose (n-m+r-k-1) (r-k) := sorry

theorem test_027_1(n:ℕ)(hn:n≥ 1):∑ k ∈ Finset.range (n+1),k*(Nat.choose n k)^2=n*Nat.choose (2*n-1) (n-1) := sorry

theorem test_028_1(n:ℕ)(hn:n≥ 1):∑ k ∈ Finset.Ico 1 (n+1),1/k*Nat.choose (2*(k-1)) (k-1) *Nat.choose (2*(n-k)) (n-k)=1/2*Nat.choose (2*n) n := sorry

theorem test_029_1 (n m:ℕ)(hn:n≥ 1) (h:n≥m): ∑ k ∈ Finset.range (m+1),(-1:ℝ)^k * Nat.choose n k = (-1:ℝ)^m * Nat.choose (n-1) m := sorry

theorem test_030_1(n:ℕ):∑ k ∈ Finset.range (n + 1),2^(n-k)*Nat.choose (n+k) (2*k)=(2^(2*n+1)+1)/3 := sorry

theorem test_031_1 (n m :ℕ) : ∑ k ∈ Finset.range (n+1), Nat.choose m k  * Nat.choose (m+n-1) k = (m+n)/n := sorry

theorem test_032_1(n:ℕ):∑'k :ℕ,(-1:ℝ)^k*Nat.choose (n-k) (m-k)*Nat.choose n k=0 := sorry

theorem test_033_1(n:ℕ):∑'k :ℕ,Nat.choose (n-k) (n-m)*Nat.choose n k=2^m*Nat.choose n m := sorry

theorem test_034_1(n:ℕ)(x:ℝ):∑ i in Finset.range (n+1),(Nat.choose n i)^2*x^i=∑ i in Finset.range (n+1),Nat.choose n i *Nat.choose (2*n-i) n*(x-1)^i := sorry

theorem test_035_1(n k:ℕ)(hn:n≥1)(hr:r≥1):∑ r in Finset.range (n+2),Nat.choose (n+k-r) (k-1)=Nat.choose (n+k+1) k := sorry

theorem test_036_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * Nat.choose n k * 2 ^ (2 * k) * (Nat.choose (2 * k) k : ℝ)⁻¹ =
      1 / (1 - 2 * n) := sorry

theorem test_037_1 (n : ℕ) (hn : n ≥ 4) :
  ∑ k ∈ Finset.range (n + 1), k ^ 4 * Nat.choose n k =
    Nat.choose n 1 * 2 ^ (n - 1) + 14 * Nat.choose n 2 * 2 ^ (n - 2) +
    36 * Nat.choose n 3 * 2 ^ (n - 3) + 24 * Nat.choose n 4 * 2 ^ (n - 4) := sorry

theorem test_038_1 (n : ℕ) (hn : 1 ≤ n):
    ∑ k ∈ Finset.range (n + 1), (4 * n).choose (4 * k) =
      (1 / 4 : ℝ) * (2 ^ (4 * n) + (-1) ^ n * 2 ^ (2 * n + 1)) := sorry

theorem test_039_1 (n m : ℕ) :
    ∑ k ∈ Finset.range (n + 1), choose (m - k) (n - k) = choose (m + 1) n := sorry

theorem test_040_1 (n : ℕ) (hn : 1 ≤ n):
  ∑ k ∈ Finset.range (n + 1), ∑ l in Finset.range (k + 1), (-1 : ℝ) ^ l * Nat.choose k l * (l + 1) ^ n = 0 := sorry

theorem test_041_1 (n : ℕ) (x : ℝ):
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k / (k + 1) * Nat.choose n k * ((1 + x) ^ (k + 1) - 1) =
    (-1) ^ n * x ^ (n + 1) / (n + 1) := sorry

theorem test_042_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n / 2 + 1), ((n + 1 - 2 * k) ^ 2) / (n + 1 - k) * Nat.choose n k = 2 ^ n := sorry

theorem test_043_1 (n : ℕ) :
  ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * choose n k * (∑ l in range k, (1 / l ^ 2 : ℝ))
  = -1 / n * ∑ k ∈ Finset.range (n + 1), 1 / k := sorry

theorem test_044_1 (n : ℕ) :
    ∑ k ∈ Finset.range n, (-1 : ℝ) ^ k * (cos (k * π / n)) ^ n =
      n / 2 ^ (n - 1) := sorry

theorem test_045_1 (n k : ℕ)(hn : n > 0)(hk : k > 0) :
    (Nat.choose n k : ℝ) = (n : ℝ) / (k : ℝ) * Nat.choose (n - 1) (k - 1) := sorry

theorem test_046_1 (n m : ℕ)(h : m ≥ 2 * n) (hn : 1 ≤ n):
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * Nat.choose n k * Nat.choose (m - 2 * k) (n - 1) = 0 := sorry

theorem test_047_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n + 1), Nat.choose (2 * n) k = 2 ^ (n - 1) + (Nat.choose (2 * n) n) / 2 := sorry

theorem test_048_1 (n m k : ℕ)(hmk : m ≤ k) (hkn : k ≤ n) :
    Nat.choose n k * Nat.choose k m = Nat.choose n (k - m) * Nat.choose (n - k + m) m := sorry

theorem test_049_1 (n : ℕ) (x : ℝ) (hn : 1 ≤ n):
    ∑ k ∈ Finset.range (n / 2 + 1), (Nat.choose n (2 * k)) * sin (2 * k * x) =
      2 ^ (n - 1) * (sin (n * x / 2) * (cos (x / 2)) ^ n + sin ((x + π) / 2) * (cos ((x + π) / 2)) ^ n) := sorry

theorem test_050_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n / 2 + 1), Nat.choose (n - k) k * 6 ^ k =
    (1 / 5) * (3 ^ (n + 1) + (-1 : ℝ) ^ n * 2 ^ (n + 1)) := sorry

theorem test_051_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n / 2 + 1), Nat.choose (n - k) k * 2 ^ k =
    (1 / 3 : ℝ) * (2 ^(n + 1) + (-1 : ℝ) ^ n) := sorry

theorem test_052_1 (n m : ℕ) (hmn : 0 < m ∨ 0 < n) :
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * (Nat.choose n k) * (m / (m + k)) =
    1 / (Nat.choose (m + n) n : ℝ) := sorry

theorem test_053_1 (n : ℕ) (x : ℝ) :
  ∑ k ∈ Finset.range (n + 1), (-1) ^ k * choose n k * (x + n - k) ^ n
  = n ! := sorry

theorem test_054_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * 2 ^ (2 * n - 2 * k) * Nat.choose (2 * n - k) k = 2 * n + 1 := sorry

theorem test_055_1 (n : ℕ) (x : ℝ) :
  ∑ k ∈ Finset.range (n + 1), (-1) ^ k * choose n k * (x - k) ^ (n + 1)
  = (x - n / 2) * (n + 1) ! := sorry

theorem test_056_1 (n m : ℕ) (hn : 1 ≤ n):
    ∑ k ∈ Finset.range (m + 1), Nat.choose (k + n - 1) (n - 1) = Nat.choose (n + m) n := sorry

theorem test_057_1 (n : ℕ) (hn : n ≥ 2) :
    ∑ k ∈ Finset.range (n - 1), Nat.choose n k = 2 ^ (n - 1) := sorry

theorem test_058_1 (n : ℕ) :
    ∑ k ∈ Ico 1 (n + 1), ((-1 : ℝ) ^ (k - 1) * (2 ^ (2 * k : ℝ) / (2 * k - 1 : ℝ)) * (Nat.choose n k : ℝ) / (Nat.choose (2 * k - 2) (k - 1)  : ℝ))=
     4 * ∑ k ∈ Finset.range n, (1 / (2 * k + 1)) := sorry

theorem test_059_1 (n : ℕ) :
    ∑ k ∈  Finset.range (n + 1), Nat.choose (n + k) (2 * k) * 2 ^ (n - k) =
    (1 / 3 : ℝ) * (2 ^ (2 * n + 1) + 1) := sorry

theorem test_060_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n + 1),
      (-1 : ℝ) ^ k / (2 * k + 1) * Nat.choose n k * (Nat.choose (2 * k) k : ℝ)⁻¹ * 2 ^ (2 * k) =
      1 / (2 * n + 1) := sorry

theorem test_061_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * 2 ^ (2 * n - 2 * k) * Nat.choose (2 * n - k + 1) k = n + 1 := sorry

theorem test_062_1 (n : ℕ) (x : ℝ):
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * Nat.choose n k * (∑ k ∈ Finset.range (n + 1), (x ^ k) / k) =
    ((-1 : ℝ) / n) * (1 - (1 - x) ^ n) := sorry

theorem test_063_1 (n p: ℕ) (x : ℝ) :
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * Nat.choose (n + p) (k + p) =
     Nat.choose (n + p - 1) n := sorry

theorem test_064_1 (n : ℕ) (hn : 1 ≤ n):
    ∑ k ∈ Finset.range (n / 2 + 1),
      (-1 : ℝ) ^ k * (n / (n - k)) * (Nat.choose (n - k) k) * 2 ^ (n - 2 * k) = 2 := sorry

theorem test_065_1 (n m r : ℕ) :
    ∑ k ∈ Finset.range (n + 1), choose n k * choose m (k + r) = choose (n + m) (n + r) := sorry

theorem test_066_1 (n m : ℕ) :
    ∑ k ∈ Finset.range (n + 1), choose n k * choose m k * k = n * choose (n + m - 1) n := sorry

theorem test_067_1 (n m r : ℕ) :
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * choose n k * choose (m + k) (r + k) = (-1 : ℝ) ^ n * choose m (n + r) := sorry

theorem test_068_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * choose n k * 10 ^ k = (-9 : ℝ) ^ n := sorry

theorem test_069_1 (n : ℕ) (hn : 2 < n) :
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * k ^ 2 * choose n k = 0 := sorry

theorem test_070_1 (n : ℕ) (hn : 1 ≤ n):
    ∑ k ∈ Finset.range (n + 1), choose (2 * n) k = 2 ^ (2 * n - 1) + (1 / 2 : ℝ) * choose (2 * n) n := sorry

theorem test_071_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n + 1), ((k : ℝ) + 1) ^ 2 * choose n k = (2 : ℝ) ^ ((n : ℝ) - 2) * (n ^ 2 + 5 * n + 4) := sorry

theorem test_072_1 (n m s : ℕ) :
    ∑ k ∈ Finset.range s, choose k n * choose (k + m) m = choose s n * choose (s + m) m * (s - n) / (m + n + 1 : ℝ) := sorry

theorem test_073_1 (x n : ℕ) :
    ∑ k ∈ Finset.range (n + 1),(Nat.choose n k) ^ 2 * choose (x + k) (2 * n)=(Nat.choose x n) ^ 2 := sorry

theorem test_074_1 (x:ℝ)(n:ℕ):(x-1)^n=∑ k ∈ Finset.range (n+1),Nat.choose n k *x^k*(-1:ℝ )^(n-k) := sorry

theorem test_075_1 (n r : ℕ) : r ! ∣ ((n+1).ascFactorial r) := sorry

theorem test_076_1 (n m k : ℕ) :
    ∑ r ∈ Finset.range (k + 1), choose n r * choose m (k - r) = choose (n + m) k := sorry

theorem test_077_1(n:ℕ): ∑ k ∈ Finset.Ico 1 (2*n),(-1:ℝ)^(k-1)*1/(Nat.choose (2*n) k)=1/(n+1) := sorry

theorem test_078_1 : catalan 5 = 42 := sorry

theorem test_079_1 (n m:ℕ):∑ k ∈ Finset.range (m+1),Nat.choose (n+k+1) k =Nat.choose (n+m+2) m := sorry

theorem test_080_1 (x : ℝ) (k : ℕ) :
  ∑' n : Set.Ici k, deriv (fun (x:ℝ) => (choose n k) * (x ^ (n : ℕ))) x
  =  ∑' n : Set.Ici k, (choose n k) * n * (x ^ ((n-1) : ℕ)) := sorry

theorem test_081_1 (n:ℕ):∑ k ∈ Finset.range (n+2),k=Nat.choose (n+2) 2 := sorry

theorem test_082_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n + 1), 2 ^ k * Nat.choose n k = 3 ^ n := sorry

theorem test_083_1 (n k : ℕ) : choose (n + 2) (k + 3) = choose (n + 1) (k + 2) + choose (n + 1) (k + 3) := sorry

theorem test_084_1 (n k:ℕ)(hn:n≥3)(hk:k≥1):
  Nat.choose n k - Nat.choose (n-3) k=Nat.choose (n-1) (k-1) +Nat.choose (n-2) (k-1)+Nat.choose (n-3) (k-1) := sorry

theorem test_085_1 (n k : ℕ) :
    n.ascFactorial (k + 1) = (k + 1) ! * (n + k).choose (k + 1) := sorry

theorem test_086_1 (n : ℕ) :
    (n + 2).descFactorial 2 = (n + 1).ascFactorial 2 := sorry

theorem test_087_1 (n : ℕ) : 2 * (2 * n)! ≤ 2 * (2 * n) ^ n * n ! := sorry

theorem test_088_1 {n k s : ℕ} (hkn : k ≤ n) (hsk : s + 1 ≤ k) :
    n.choose k * k.choose (s + 1) = n.choose (s + 1) * (n - s - 1).choose (k - s - 1) := sorry

theorem test_089_1 (n m r : ℕ) :
    ∑ k ∈ Finset.range (n + 1), choose m k * choose r (n - k) = choose (m + r) n := sorry

theorem test_090_1(n:ℕ) : 3^n = ∑ k ∈ Finset.range (n+1), Nat.choose n k *2^k := sorry

theorem test_091_1(n:ℕ) : 2^n = ∑ k ∈ Finset.range (n+1),(-1:ℝ)^k*Nat.choose n k *3^(n-k) := sorry

theorem test_092_1 (n : ℕ) : choose n (2 * n + 2) = 0 := sorry

theorem test_093_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n / 2 + 1), (-1 : ℝ) ^ k * Nat.choose (n - k ) k =(2 / Real.sqrt 3) *  Real.sin ((n + 1) * Real.pi / 3) := sorry

theorem test_094_1 (i j : ℕ) :
    (i + j + 1).choose (j + 1) * i ! * (j + 1) ! = (i + (j + 1))! := sorry

theorem test_095_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n + 1), (-1 : ℝ) ^ k * 1/ (Nat.choose n k) =
      (n + 1) / (n + 2) * (1 + (-1) ^ n) := sorry

theorem test_096_1 (n : ℕ) (hn : 3 ≤ n): Nat.choose n 3 = n ! / (3! * (n - 3)!) := sorry

theorem test_097_1(n:ℕ):∑ k ∈ Finset.range (n + 1),(Nat.fib k)^2=Nat.fib n * Nat.fib (n+1) := sorry

theorem test_098_1 (n : ℕ) :
    ∑ k ∈ Finset.range (n / 2 + 1), Nat.choose (n - k) k =
    (1 / Real.sqrt 5) * ( ((1 + Real.sqrt 5) / 2) ^ (n + 1) - ((1 - Real.sqrt 5) / 2) ^ (n + 1)) := sorry

theorem test_099_1 (n k : ℕ) (hk : 1 < k) :
    choose (n + 1) (k - 1) = choose n (k - 2) + choose n (k - 1) := sorry

theorem test_100_1 (n:ℕ) : ∑ m ∈ range (n + 1), (-1:ℝ) ^ m * 10 ^ (n- m) * n.choose m=9^n := sorry
