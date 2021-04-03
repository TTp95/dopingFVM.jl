"""

"""
function check_folder(
    path::String,
)

    if isdir(path)
        print("Folder \"$path\" is already created... \n")
    else
        mkdir(path)
        print("Folder \"$path\" has been created...\n")
    end

    return nothing
end
