import Mathlib

open Nat Finset BigOperators

def stirling_first : ℕ → ℕ → ℕ
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | _ + 1, 0 => 0
  | n + 1, k + 1 => n * stirling_first n (k + 1) + stirling_first n k
#eval stirling_first 4 2


theorem stirling_first_zero_right : stirling_first 0 0 = 1 := by simp [stirling_first]

theorem stirling_first_zero_succ (k : ℕ) : stirling_first 0 (k + 1) = 0 := by simp [stirling_first]

theorem stirling_first_succ_zero (n : ℕ) : stirling_first (n + 1) 0 = 0 := by simp [stirling_first]

theorem stirling_first_succ_left (n k : ℕ) (hk : 0 < k) :
    stirling_first (n + 1) k = n * stirling_first n k + stirling_first n (k - 1) := by
  obtain ⟨l, rfl⟩ : ∃ l, k = l + 1 := Nat.exists_eq_add_of_le' hk
  rfl

theorem stirling_first_succ_right (n k : ℕ) (hn : 0 < n) :
    stirling_first n (k + 1) = (n - 1) * stirling_first (n - 1) (k + 1) + stirling_first (n - 1) k := by
  obtain ⟨l, rfl⟩ : ∃ l, n = l + 1 := Nat.exists_eq_add_of_le' hn
  rfl


theorem stirling_first_succ_succ (n k : ℕ) :
    stirling_first (n + 1) (k + 1) = n * stirling_first n (k + 1) +  stirling_first n k := by
  rw [stirling_first]

theorem stirling_first_eq_zero_of_lt : ∀ {n k :ℕ }, n < k → stirling_first n k= 0
  | _, 0, hk => absurd hk (Nat.not_lt_zero _)
  | 0, _ + 1, _ => by rw [stirling_first ]
  | n + 1, k + 1, hk => by
    have hnk : n < k := lt_of_succ_lt_succ hk
    have hnk1 : n < k + 1 := lt_of_succ_lt hk
    rw [stirling_first_succ_succ, stirling_first_eq_zero_of_lt hnk ,stirling_first_eq_zero_of_lt hnk1]
    rfl

theorem stirling_first_self (n : ℕ) : stirling_first n n = 1 := by
  induction n <;> simp [*, stirling_first, stirling_first_eq_zero_of_lt (lt_succ_self _)]

theorem stirling_first_succ_self_left (n : ℕ) : stirling_first (n + 1) n = (n * (n + 1)) / 2 := by
  induction' n with n hn
  · simp [stirling_first]
  · have h₀ : (n + 1) = (2 * (n + 1)) / 2 := by
      rw [mul_comm, Nat.mul_div_assoc, Nat.div_self, mul_one]
      · omega
      · exact Nat.dvd_refl _
    rw [stirling_first_succ_succ, hn, stirling_first_self, mul_one]
    nth_rw 1 [h₀]
    rw [← Nat.add_div_of_dvd_left]
    · ring_nf
    · suffices h₁ : Even (n * (n + 1)) from by
        rw [even_iff_two_dvd] at h₁
        exact h₁
      exact Nat.even_mul_succ_self n


theorem stirling_first_one_right (n : ℕ) : stirling_first (n + 1) 1 = n ! := by
  induction' n with n hn
  · simp [stirling_first]
  · rw [stirling_first_succ_succ, zero_add, hn, stirling_first_succ_zero]
    simp [Nat.sub_self, Nat.factorial_succ]

def stirling_second : ℕ → ℕ → ℕ
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | _ + 1, 0 => 0
  | n + 1, k + 1 =>
     (k + 1) * stirling_second n (k + 1) + stirling_second n k



theorem stirling_second_zero_right (n : ℕ) : stirling_second 0 0 = 1 := by cases n <;> simp [stirling_second]

theorem stirling_second_zero_succ (k : ℕ) : stirling_second 0 (k + 1) = 0 := by simp [stirling_second]

theorem stirling_second_succ_left (n k : ℕ) (hk : 0 < k) :
    stirling_second (n + 1) k = k * stirling_second n k + stirling_second n (k - 1) := by
  obtain ⟨l, rfl⟩ : ∃ l, k = l + 1 := Nat.exists_eq_add_of_le' hk
  rfl

theorem stirling_second_succ_right (n k : ℕ) (hn : 0 < n) :
    stirling_second n (k + 1) = (k + 1) * stirling_second (n - 1) (k + 1) + stirling_second (n - 1) k := by
  obtain ⟨l, rfl⟩ : ∃ l, n = l + 1 := Nat.exists_eq_add_of_le' hn
  rfl


theorem stirling_second_succ_succ (n k : ℕ) :
    stirling_second (n + 1) (k + 1) =  (k + 1) * stirling_second n (k + 1) + stirling_second n k := by
  rw [stirling_second]

theorem stirling_second_eq_zero_of_lt : ∀ {n k :ℕ }, n < k → stirling_second n k = 0
  | _, 0, hk => absurd hk (Nat.not_lt_zero _)
  | 0, _ + 1, _ => by rw [stirling_second ]
  | n + 1, k + 1, hk => by
    have hnk : n < k := lt_of_succ_lt_succ hk
    have hnk1 : n < k + 1 := lt_of_succ_lt hk
    rw [stirling_second_succ_succ, stirling_second_eq_zero_of_lt hnk ,stirling_second_eq_zero_of_lt hnk1]
    simp

theorem stirling_second_zero_right_n (n : ℕ) (hn : 0 < n): stirling_second n 0 = 0 := by
  have h : 1 ≤ n := Nat.succ_le_of_lt hn
  rw [show n = n - 1 + 1 by rw [Nat.sub_add_cancel h]]
  simp [stirling_second]


theorem stirling_second_zero_right' (n : ℕ) : stirling_second (n+1) 0 = 0 := by cases n <;> simp [stirling_second]

theorem stirling_second_one_right (n : ℕ) : stirling_second (n+1) 1 = 1 := by
  simp [stirling_second]
  induction' n with n ih
  . simp [stirling_second]
  . rw [stirling_second_zero_right']
    simp
    nth_rw 2 [ show 1 = 0 + 1 by ring ]
    rw [stirling_second_succ_succ]
    simp
    exact ih
