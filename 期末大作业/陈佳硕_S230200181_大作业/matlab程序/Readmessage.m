%%%��ȡinp�ļ�ģ����Ϣ����
function [Nodes, Elements] = Readmessage( fname )
fid = fopen(fname,'rt');  %fname�ļ���   r��ȡ  t��txt��ʽ��
S = textscan(fid,'%s','Delimiter','\n'); %�Ѿ��򿪵��ļ�  �ַ�������ʽ��ȡ   'Delimiter','\n'�ָ���Ϊ���з�  Ĭ�Ϸָ���Ϊ�ո�
S = S{1};
%�ҵ�Node�ؼ������ڵ�λ��
idxS = strfind(S, 'Node');  %����Ԫ������ ��������û����ӦԪ�أ��򷵻ؿ�
idx1 = find(not(cellfun(@isempty, idxS))); %cellfun(fun,A) ��Ԫ������A�ֱ�ʹ�ú���fun   isempty(A) AΪ�շ����߼�ֵ1   find Ѱ�ҷ�0Ԫ�ص�����
%�ҵ�Element�ؼ������ڵ�λ��
idxS = strfind(S, 'Element');
idx2 = find(not(cellfun(@isempty, idxS)));
%�ҵ�Nset�ؼ�������λ��
idxS = strfind(S, 'Nset');
idx3 = find(not(cellfun(@isempty, idxS)));
% ȡ���ڵ���Ϣ(Ԫ������)
Nodes = S(idx1(1)+1:idx2(1)-1);  %��Ԫ��������ʽȡ��
%��Ԫ������ת��Ϊ����
Nodes = cell2mat(cellfun(@str2num,Nodes,'UniformOutput',false));  %'UniformOutput',false  ��Ԫ����ʽ�������ֵ
% ȡ����Ԫ(Ԫ������)
elements = S(idx2+1:idx3(1)-1) ;
% ��Ԫ������ת��Ϊ����
Elements = cell2mat(cellfun(@str2num,elements,'UniformOutput',false));
Nodes=Nodes(:,2:end);
Elements=Elements(:,2:end);
end