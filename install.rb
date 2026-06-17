#!/usr/bin/env ruby
# frozen_string_literal: true

# FLIT installer.
#
# Common usage:
#   ruby install.rb
#   ruby install.rb --without-hdf5
#   ruby install.rb --intel-root ~/intel --hdf5-root ~/tool/hdf5 --hdf5-version 2.1.1 -j 16
#
# Environment overrides:
#   FLIT_USE_HDF5=true|false
#   FLIT_INTEL_ROOT=/path/to/intel
#   FLIT_HDF5_ROOT=/path/to/hdf5
#   FLIT_HDF5_VERSION=2.1.1
#   FLIT_HDF5_URL=https://.../hdf5-2.1.1.tar.gz
#   FLIT_HDF5_SOURCE=/path/to/hdf5-source-or-archive
#   FLIT_JOBS=16
#   FLIT_CLEAN=auto|true|false
#   SHELL_CONFIG=/path/to/shell/startup/file
#   FC=ifx CC=icx CXX=icpx MAKE=make CMAKE=cmake

require "etc"
require "fileutils"
require "open-uri"
require "optparse"
require "shellwords"
require "uri"

REPO_ROOT = File.expand_path(__dir__)
SRC_DIR = File.join(REPO_ROOT, "src")
MAKEFILE_IN = File.join(SRC_DIR, "Makefile.in")

def env_bool(name, default)
  value = ENV[name]
  return default if value.nil? || value.strip.empty?

  case value.strip.downcase
  when "1", "yes", "y", "true", "on"
    true
  when "0", "no", "n", "false", "off"
    false
  else
    raise ArgumentError, "#{name} must be true/false, on/off, yes/no, or 1/0"
  end
end

def clean_policy(value)
  case value.to_s.strip.downcase
  when "", "auto"
    :auto
  when "1", "yes", "y", "true", "on", "always"
    :always
  when "0", "no", "n", "false", "off", "never"
    :never
  else
    raise ArgumentError, "FLIT_CLEAN must be auto, true, false, always, or never"
  end
end

def expand_path(path)
  value = path.to_s
  home = ENV["HOME"]

  if home && value.match?(%r{\A(~|\$HOME|\$\{HOME\}|\$\(HOME\))(/|\z)})
    value = value.sub(%r{\A(~|\$HOME|\$\{HOME\}|\$\(HOME\))}, home)
  end

  File.expand_path(value)
end

def make_path(path)
  expanded = expand_path(path)
  home = ENV["HOME"] && File.expand_path(ENV["HOME"])
  return expanded if home.nil? || home.empty?

  if expanded == home
    "$(HOME)"
  elsif expanded.start_with?(home + File::SEPARATOR)
    "$(HOME)/#{expanded.delete_prefix(home + File::SEPARATOR)}"
  else
    expanded
  end
end

def default_hdf5_url(version)
  #"https://support.hdfgroup.org/releases/hdf5/#{version}/downloads/hdf5-#{version}.tar.gz"
  "https://github.com/HDFGroup/hdf5/releases/download/#{version}/hdf5-#{version}.tar.gz"
end

def executable_available?(command)
  return File.executable?(command) if command.include?(File::SEPARATOR)

  ENV.fetch("PATH", "").split(File::PATH_SEPARATOR).any? do |dir|
    File.executable?(File.join(dir, command))
  end
end

def run_command(*command, chdir: nil)
  location = chdir ? " in #{chdir}" : ""
  puts "==> #{command.shelljoin}#{location}"

  env = chdir ? { "PWD" => chdir } : {}

  options = {}
  options[:chdir] = chdir if chdir

  ok = system(env, *command, **options)
  raise "Command failed: #{command.shelljoin}" unless ok
end

def replace_make_setting(text, key, value)
  pattern = /^(\s*#{Regexp.escape(key)}\s*(?:\?=|=)\s*).*$/
  raise "Could not find '#{key}' in #{MAKEFILE_IN}" unless text.match?(pattern)

  text.sub(pattern) { "#{$1}#{value}" }
end

def configure_flit_makefile(use_hdf5:, intel_root:, hdf5_root:)
  old_text = File.read(MAKEFILE_IN)
  text = old_text
  text = replace_make_setting(text, "use_hdf5", use_hdf5 ? "on" : "off")
  text = replace_make_setting(text, "intelroot", make_path(intel_root))
  text = replace_make_setting(text, "hdf5root", make_path(hdf5_root))
  text = text.gsub(%r{-L/?\$\((?:HOME|hdf5root)\)/(?:tool/)?hdf5/lib}, "-L$(hdf5root)/lib")

  changed = text != old_text
  File.write(MAKEFILE_IN, text) if changed
  changed
end

def hdf5_installed?(root)
  root = expand_path(root)
  include_dir = File.join(root, "include")
  lib_dirs = [File.join(root, "lib"), File.join(root, "lib64")]

  has_fortran_module = !Dir[File.join(include_dir, "**", "hdf5.mod")].empty?
  has_fortran_library = lib_dirs.any? do |dir|
    !Dir[File.join(dir, "libhdf5_fortran.{a,so,dylib}")].empty?
  end
  has_c_library = lib_dirs.any? do |dir|
    !Dir[File.join(dir, "libhdf5.{a,so,dylib}")].empty?
  end
  has_wrapper = %w[h5fc h5pfc].any? do |program|
    File.executable?(File.join(root, "bin", program))
  end

  has_wrapper || (has_fortran_module && has_fortran_library && has_c_library)
end

def prepend_env_path(name, path)
  current = ENV[name].to_s
  paths = current.empty? ? [] : current.split(File::PATH_SEPARATOR)
  paths.unshift(path) unless paths.include?(path)
  ENV[name] = paths.join(File::PATH_SEPARATOR)
end

def apply_hdf5_env(root)
  ENV["HDF5_ROOT"] = root
  prepend_env_path("PATH", File.join(root, "bin"))
  prepend_env_path("CPATH", File.join(root, "include"))
  prepend_env_path("LD_LIBRARY_PATH", File.join(root, "lib"))
end

def append_hdf5_shell_setup(root, shell_config)
  shell = File.basename(ENV["SHELL"].to_s)
  block_begin = "# >>> FLIT HDF5 paths >>>"
  block_end = "# <<< FLIT HDF5 paths <<<"

  if File.file?(shell_config) && File.read(shell_config).include?(block_begin)
    puts "HDF5 environment setup already exists in #{shell_config}"
    return
  end

  snippet =
    case shell
    when "csh", "tcsh"
      <<~CSH

        #{block_begin}
        setenv HDF5_ROOT "#{root}"
        setenv PATH "#{File.join(root, "bin")}:$PATH"
        if ( $?CPATH ) then
            setenv CPATH "#{File.join(root, "include")}:$CPATH"
        else
            setenv CPATH "#{File.join(root, "include")}"
        endif
        if ( $?LD_LIBRARY_PATH ) then
            setenv LD_LIBRARY_PATH "#{File.join(root, "lib")}:$LD_LIBRARY_PATH"
        else
            setenv LD_LIBRARY_PATH "#{File.join(root, "lib")}"
        endif
        #{block_end}
      CSH
    else
      escaped_root = Shellwords.escape(root)
      escaped_bin = Shellwords.escape(File.join(root, "bin"))
      escaped_include = Shellwords.escape(File.join(root, "include"))
      escaped_lib = Shellwords.escape(File.join(root, "lib"))

      <<~SH

        #{block_begin}
        export HDF5_ROOT=#{escaped_root}
        export PATH=#{escaped_bin}:$PATH
        export CPATH=#{escaped_include}:${CPATH:-}
        export LD_LIBRARY_PATH=#{escaped_lib}:${LD_LIBRARY_PATH:-}
        #{block_end}
      SH
    end

  FileUtils.mkdir_p(File.dirname(shell_config))
  File.open(shell_config, "a") { |file| file.write(snippet) }

  puts "Added HDF5 environment setup to #{shell_config}"
  puts "Restart your shell or run: source #{shell_config}"
end

def ensure_required_tools!(*commands)
  missing = commands.compact.reject { |command| executable_available?(command) }
  return if missing.empty?

  raise <<~MSG
    Missing required command(s): #{missing.join(", ")}
    Load your compiler environment first, for example:
      source ~/tool/intel/oneapi/setvars.sh
    Then run install.rb again.
  MSG
end

def safe_remove_dir(path, allowed_parent)
  path = expand_path(path)
  allowed_parent = expand_path(allowed_parent)

  unless path.start_with?(allowed_parent + File::SEPARATOR)
    raise "Refusing to remove #{path}; it is outside #{allowed_parent}"
  end

  FileUtils.rm_rf(path)
end

def download_hdf5_source(options)
  source = options[:hdf5_source]
  return expand_path(source) if source && !source.strip.empty?

  url = options[:hdf5_url] || default_hdf5_url(options[:hdf5_version])
  uri = URI.parse(url)
  filename = File.basename(uri.path)
  raise "Cannot determine archive filename from #{url}" if filename.empty?

  archive = File.join(options[:build_root], filename)
  return archive if File.file?(archive) && File.size(archive).positive?

  FileUtils.mkdir_p(options[:build_root])
  puts "==> Downloading #{url}"
  URI.open(uri, "rb") do |remote|
    File.open(archive, "wb") do |file|
      IO.copy_stream(remote, file)
    end
  end

  archive
end

def prepare_hdf5_source(source, options)
  return expand_path(source) if File.directory?(source)

  ensure_required_tools!("tar")

  build_root = expand_path(options[:build_root])
  version = options[:hdf5_version]

  extract_dir = File.join(build_root, "hdf5-#{version}")

  safe_remove_dir(extract_dir, build_root)
  FileUtils.mkdir_p(extract_dir)

  archive = expand_path(source)

  puts "Extracting #{archive}"
  puts "Into #{extract_dir}"

  run_command(
    "tar", "-xf", archive,
    "--strip-components=2",
    "-C", extract_dir
  )

  extract_dir
end

def install_hdf5(options)
  ensure_required_tools!(options[:cmake], options[:fc], options[:cc], options[:cxx])

  source = download_hdf5_source(options)
  source_dir = prepare_hdf5_source(source, options)
  build_dir = File.join(options[:build_root], "hdf5-#{options[:hdf5_version]}-build")

  safe_remove_dir(build_dir, options[:build_root])
  FileUtils.mkdir_p(build_dir)
  FileUtils.mkdir_p(options[:hdf5_root])

  cmake_args = [
    options[:cmake],
    "-S", source_dir,
    "-B", build_dir,
    "-DCMAKE_Fortran_COMPILER=#{options[:fc]}",
    "-DCMAKE_C_COMPILER=#{options[:cc]}",
    "-DCMAKE_CXX_COMPILER=#{options[:cxx]}",
    "-DCMAKE_Fortran_FLAGS=-O2 -w",
    "-DCMAKE_C_FLAGS=-O2 -w",
    "-DCMAKE_CXX_FLAGS=-O2 -w",
    "-DCMAKE_BUILD_TYPE=#{options[:cmake_build_type]}",
    "-DCMAKE_INSTALL_PREFIX=#{options[:hdf5_root]}",
    "-DCMAKE_INSTALL_LIBDIR=lib",
    "-DHDF5_BUILD_FORTRAN=ON",
    "-DHDF5_BUILD_EXAMPLES=OFF",
    "-DBUILD_TESTING=OFF"
  ]

  run_command(*cmake_args)
  run_command(options[:cmake], "--build", build_dir, "-j#{options[:jobs]}")
  run_command(options[:cmake], "--install", build_dir)

  unless hdf5_installed?(options[:hdf5_root])
    raise "HDF5 install completed, but FLIT could not find the Fortran module/library in #{options[:hdf5_root]}"
  end

end

def parse_options
  options = {
    use_hdf5: env_bool("FLIT_USE_HDF5", true),
    intel_root: expand_path(ENV.fetch("FLIT_INTEL_ROOT", "~/intel")),
    hdf5_root: expand_path(ENV.fetch("FLIT_HDF5_ROOT", "~/tool/hdf5")),
    hdf5_version: ENV.fetch("FLIT_HDF5_VERSION", "2.1.1"),
    hdf5_url: ENV["FLIT_HDF5_URL"],
    hdf5_source: ENV["FLIT_HDF5_SOURCE"],
    build_root: expand_path(ENV.fetch("FLIT_BUILD_ROOT", File.join(REPO_ROOT, ".flit-build"))),
    jobs: Integer(ENV.fetch("FLIT_JOBS", Etc.nprocessors.to_s)),
    cmake: ENV.fetch("CMAKE", "cmake"),
    make: ENV.fetch("MAKE", "make"),
    fc: ENV.fetch("FC", "ifx"),
    cc: ENV.fetch("CC", "icx"),
    cxx: ENV.fetch("CXX", "icpx"),
    cmake_build_type: ENV.fetch("FLIT_CMAKE_BUILD_TYPE", "Release"),
    clean_policy: clean_policy(ENV.fetch("FLIT_CLEAN", "auto")),
    skip_flit_build: env_bool("FLIT_SKIP_BUILD", false),
    shell_config: expand_path(ENV.fetch("SHELL_CONFIG", "~/.bashrc"))
  }

  OptionParser.new do |parser|
    parser.banner = "Usage: ruby install.rb [options]"

    parser.on("--with-hdf5", "Build FLIT with HDF5, installing HDF5 if needed") do
      options[:use_hdf5] = true
    end

    parser.on("--without-hdf5", "Build FLIT without HDF5") do
      options[:use_hdf5] = false
    end

    parser.on("--intel-root PATH", "Intel compiler/MKL/MPI root used by src/Makefile.in") do |value|
      options[:intel_root] = expand_path(value)
    end

    parser.on("--hdf5-root PATH", "HDF5 installation prefix") do |value|
      options[:hdf5_root] = expand_path(value)
    end

    parser.on("--hdf5-version VERSION", "HDF5 source version to download") do |value|
      options[:hdf5_version] = value
    end

    parser.on("--hdf5-url URL", "Download URL for the HDF5 source archive") do |value|
      options[:hdf5_url] = value
    end

    parser.on("--hdf5-source PATH", "Existing HDF5 source directory or archive") do |value|
      options[:hdf5_source] = expand_path(value)
    end

    parser.on("--build-root PATH", "Temporary build/download directory") do |value|
      options[:build_root] = expand_path(value)
    end

    parser.on("--shell-config PATH", "Shell startup file to patch with HDF5 paths") do |value|
      options[:shell_config] = expand_path(value)
    end

    parser.on("-j", "--jobs N", Integer, "Parallel build jobs") do |value|
      options[:jobs] = value
    end

    parser.on("--skip-flit-build", "Configure dependencies but do not run make in src") do
      options[:skip_flit_build] = true
    end

    parser.on("--clean", "Always run make clean before building FLIT") do
      options[:clean_policy] = :always
    end

    parser.on("--no-clean", "Do not run make clean before building FLIT") do
      options[:clean_policy] = :never
    end

    parser.on("-h", "--help", "Show this help") do
      puts parser
      exit
    end
  end.parse!

  raise "jobs must be at least 1" if options[:jobs] < 1

  options
end

def main
  options = parse_options
  FileUtils.mkdir_p(options[:build_root])

  use_hdf5_for_flit = false
  if options[:use_hdf5]
    if hdf5_installed?(options[:hdf5_root])
      puts "==> Found HDF5 in #{options[:hdf5_root]}"
    else
      puts "==> HDF5 was not found in #{options[:hdf5_root]}; installing it now"
      install_hdf5(options)
    end
    apply_hdf5_env(options[:hdf5_root])
    append_hdf5_shell_setup(options[:hdf5_root], options[:shell_config])
    use_hdf5_for_flit = true
  else
    puts "==> HDF5 disabled; FLIT will be built without HDF5 support"
  end

  makefile_changed = configure_flit_makefile(
    use_hdf5: use_hdf5_for_flit,
    intel_root: options[:intel_root],
    hdf5_root: options[:hdf5_root]
  )
  puts "==> Configured #{MAKEFILE_IN}"
  puts "    use_hdf5=#{use_hdf5_for_flit ? "on" : "off"}"
  puts "    intelroot=#{make_path(options[:intel_root])}"
  puts "    hdf5root=#{make_path(options[:hdf5_root])}"

  return if options[:skip_flit_build]

  ensure_required_tools!(options[:make])
  should_clean = case options[:clean_policy]
                 when :always
                   true
                 when :never
                   false
                 else
                   makefile_changed
                 end

  if should_clean
    reason = options[:clean_policy] == :auto ? " because #{File.basename(MAKEFILE_IN)} changed" : ""
    puts "==> Running make clean#{reason}"
    run_command(options[:make], "clean", chdir: SRC_DIR)
  else
    puts "==> Skipping make clean"
  end

  run_command(options[:make], "use_hdf5=#{use_hdf5_for_flit ? "on" : "off"}", chdir: SRC_DIR)
end

main if $PROGRAM_NAME == __FILE__
