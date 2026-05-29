#!/bin/sh

set -eu

patch_shell=false

usage() {
    cat <<'EOF'
Usage: ./set_intel_paths.sh [options]

Creates a compact Intel toolchain directory, then writes temporary
shell-specific environment files in the current directory. The sh/bash/zsh
file is sourced automatically by this script, and both temporary files are
removed before the script exits.

Environment overrides:
  FLIT_INTEL_BASE  Full Intel oneAPI installation root
                   default: $HOME/tool/intel/oneapi
  FLIT_INTEL_ROOT  Compact FLIT Intel root used by Makefile.in
                   default: $HOME/intel

Options:
  --intel-base PATH  Full Intel oneAPI installation root
  --intel-base=PATH  Same as above
  --intel-root PATH  Compact FLIT Intel root used by Makefile.in
  --intel-root=PATH  Same as above
  --patch-shell    Append the right source command to your login shell rc file
  -h, --help       Show this help
EOF
}

expand_home_path() {
    case "$1" in
        "~")
            printf '%s\n' "$HOME"
            ;;
        "~"/*)
            printf '%s\n' "$HOME/${1#~/}"
            ;;
        '$HOME')
            printf '%s\n' "$HOME"
            ;;
        '$HOME'/*)
            printf '%s\n' "$HOME/${1#\$HOME/}"
            ;;
        '${HOME}')
            printf '%s\n' "$HOME"
            ;;
        '${HOME}'/*)
            printf '%s\n' "$HOME/${1#\$\{HOME\}/}"
            ;;
        *)
            printf '%s\n' "$1"
            ;;
    esac
}

intel_base=${FLIT_INTEL_BASE:-"$HOME/tool/intel/oneapi"}
intel_root=${FLIT_INTEL_ROOT:-"$HOME/intel"}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --intel-base=*)
            intel_base=${1#*=}
            ;;
        --intel-base)
            shift
            if [ "$#" -eq 0 ]; then
                echo "--intel-base requires a path" >&2
                exit 1
            fi
            intel_base=$1
            ;;
        --intel-root=*)
            intel_root=${1#*=}
            ;;
        --intel-root)
            shift
            if [ "$#" -eq 0 ]; then
                echo "--intel-root requires a path" >&2
                exit 1
            fi
            intel_root=$1
            ;;
        --patch-shell)
            patch_shell=true
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
    shift
done

intel_base=$(expand_home_path "$intel_base")
intel_root=$(expand_home_path "$intel_root")

home_relative_path() {
    case "$1" in
        "$HOME")
            printf '%s\n' '${HOME}'
            ;;
        "$HOME"/*)
            printf '%s\n' '${HOME}'/"${1#"$HOME"/}"
            ;;
        *)
            printf '%s\n' "$1"
            ;;
    esac
}

shell_name() {
    basename "${SHELL:-sh}"
}

startup_file_for_shell() {
    case "$(shell_name)" in
        bash)
            printf '%s\n' "$HOME/.bashrc"
            ;;
        zsh)
            printf '%s\n' "$HOME/.zshrc"
            ;;
        csh)
            printf '%s\n' "$HOME/.cshrc"
            ;;
        tcsh)
            printf '%s\n' "$HOME/.tcshrc"
            ;;
        *)
            printf '%s\n' "$HOME/.profile"
            ;;
    esac
}

append_startup_block_for_shell() {
    shell_rc=$1
    block_begin="# >>> FLIT Intel toolchain paths >>>"
    block_end="# <<< FLIT Intel toolchain paths <<<"

    if grep -F "$block_begin" "$shell_rc" >/dev/null 2>&1; then
        echo "Startup file already contains FLIT Intel setup: $shell_rc"
        return
    fi

    {
        printf '\n%s\n' "$block_begin"
        case "$(shell_name)" in
            csh|tcsh)
                cat "$env_csh"
                ;;
            *)
                cat "$env_sh"
                ;;
        esac
        printf '%s\n' "$block_end"
    } >> "$shell_rc"

    echo "Added FLIT Intel setup to $shell_rc"
}

cleanup_env_files() {
    rm -f "$env_sh" "$env_csh"
}

source_temp_env() {
    . "$env_sh"
    echo "Sourced temporary environment setup: $env_sh"
}

################################################################################
# Create directories for Intel's compiler suite

rm -f "$intel_root/include" "$intel_root/lib"
rm -rf "$intel_root/bin" "$intel_root/mkl" "$intel_root/mpi"
mkdir -p "$intel_root/bin" "$intel_root/mkl" "$intel_root/mpi"

ln -sfn "$intel_base/compiler/latest/lib/" "$intel_root/lib"

ln -sfn "$intel_base/compiler/latest/bin/ifx" "$intel_root/bin/ifx"
ln -sfn "$intel_base/compiler/latest/bin/icx" "$intel_root/bin/icx"
ln -sfn "$intel_base/compiler/latest/bin/icpx" "$intel_root/bin/icpx"

ln -sfn "$intel_base/mkl/latest/include/" "$intel_root/mkl/include"
ln -sfn "$intel_base/mkl/latest/lib/intel64/" "$intel_root/mkl/lib"

ln -sfn "$intel_base/mpi/latest/bin/" "$intel_root/mpi/bin"
ln -sfn "$intel_base/mpi/latest/include/mpi" "$intel_root/mpi/include"
ln -sfn "$intel_base/mpi/latest/lib/" "$intel_root/mpi/lib"
ln -sfn "$intel_base/mpi/latest/env/" "$intel_root/mpi/env"

################################################################################
# Write shell-specific environment snippets

root_intel_rc=$(home_relative_path "$intel_root")
env_sh="$PWD/flit-intel-env.sh"
env_csh="$PWD/flit-intel-env.csh"

cat > "$env_sh" <<EOF
# Source this file from sh, bash, or zsh before building FLIT.
export root_intel=$root_intel_rc
export PATH="\${root_intel}/bin:\${root_intel}/mpi/bin:\${PATH}"
export CPATH="\${root_intel}/mkl/include:\${root_intel}/mpi/include:\${CPATH:-}"
export LD_LIBRARY_PATH="\${root_intel}/lib:\${root_intel}/mkl/lib:\${root_intel}/mpi/lib:\${root_intel}/mpi/lib/release:\${LD_LIBRARY_PATH:-}"
EOF

cat > "$env_csh" <<EOF
# Source this file from csh or tcsh before building FLIT.
setenv root_intel $root_intel_rc
setenv PATH "\${root_intel}/bin:\${root_intel}/mpi/bin:\${PATH}"
if ( \$?CPATH ) then
    setenv CPATH "\${root_intel}/mkl/include:\${root_intel}/mpi/include:\${CPATH}"
else
    setenv CPATH "\${root_intel}/mkl/include:\${root_intel}/mpi/include"
endif
if ( \$?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH "\${root_intel}/lib:\${root_intel}/mkl/lib:\${root_intel}/mpi/lib:\${root_intel}/mpi/lib/release:\${LD_LIBRARY_PATH}"
else
    setenv LD_LIBRARY_PATH "\${root_intel}/lib:\${root_intel}/mkl/lib:\${root_intel}/mpi/lib:\${root_intel}/mpi/lib/release"
endif
EOF

echo "Created Intel root: $intel_root"
echo "Wrote temporary environment files:"
echo "  $env_sh"
echo "  $env_csh"

source_temp_env

if [ "$patch_shell" = true ]; then
    shell_rc=$(startup_file_for_shell)

    touch "$shell_rc"
    append_startup_block_for_shell "$shell_rc"
else
    echo
    echo "The temporary files are only used during this setup step."
    echo "To make these paths persist in your current sh/bash/zsh session,"
    echo "source this script instead of executing it:"
    echo "  ./set_intel_paths.sh [options]"
    echo
    echo "For persistent login-shell setup, run:"
    echo "  ./set_intel_paths.sh --patch-shell [options]"
fi

cleanup_env_files
echo "Removed temporary environment files from $PWD"
