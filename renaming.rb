require("pry")



old_names = []
new_names = []

File.open("nazwy_szpiczak.txt", "r") do |f|
    f.each_line do |line|
        old_names << line.split[0]
        new_names << line.split[1]
    end
end
#binding.pry

old_names.each_with_index { |old_name, i|
    f = File.open("CEL_szpiczak/" + old_name + ".CEL")
    File.rename(f, "CEL_szpiczak/" + new_names[i] + ".CEL")
}

