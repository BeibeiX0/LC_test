import Lake
open Lake DSL

package "Theorem_Library_Comb" where
  version := v!"0.1.0"
  keywords := #["math"]
  leanOptions := #[
    ⟨`pp.unicode.fun, true⟩ -- pretty-prints `fun a ↦ b`
  ]
def leanVersion : String := s!"v{Lean.versionString}"

require mathlib from git
"https://github.com/leanprover-community/mathlib4.git" @ leanVersion


@[default_target]
lean_lib «Theorem» where
  -- add any library configuration options here
