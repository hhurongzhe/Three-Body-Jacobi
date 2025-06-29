import subprocess

# * This script is used to commit changes to a git repository.
message = "update"


command1 = "git add ."
command2 = f"git commit -m '{message}'"
command3 = "git pull origin main --rebase"
command4 = "git push origin main"
command = f"{command1} && {command2} && {command3} && {command4}"
result = subprocess.run(command, shell=True, capture_output=True, text=True)
print(result.stdout)
if result.stderr:
    print("standard error:")
    print(result.stderr)
