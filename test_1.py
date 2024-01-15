top_file_path = './topol.top'
itp_files_name = [ "inst_1.itp" , "inst_2.itp" ]



with open(top_file_path, 'r') as f:
    lines = f.readlines()

# [ molecules ] 섹션 찾기
molecules_start = lines.index('[ moleculetype ]\n')
txt_end = len(lines)

print(f"txt_end :{txt_end}")

for j in range(len(itp_files_name)):
    lines.append('#')

print(f"len(itp_files_name) :{len(itp_files_name)}")
print(len(lines))

for k in range( txt_end - molecules_start  ):
    lines[ txt_end + len(itp_files_name) - k -1 ] = lines[ txt_end - k - 1 ]

# [ molecules ] 섹션 업데이트
for i in range(len(itp_files_name)):
    lines[molecules_start - i ] ="#include " +  '\"' + itp_files_name[i] + '\"\n'
    print("#include " +  '\"' + itp_files_name[i] + '\"\n')
# top 파일 쓰기
with open(top_file_path, 'w') as f:
    f.writelines(lines)