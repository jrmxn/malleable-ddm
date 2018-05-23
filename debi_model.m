function op = debi_model(ip,ip_type,op_req)
if strcmpi(ip_type,'de')&&strcmpi(op_req,'bi')
    op = de2bi(ip,22,'left-msb');
elseif strcmpi(ip_type,'bi')&&strcmpi(op_req,'de')
    op = bi2de(ip,'left-msb');
elseif strcmpi(ip_type,'bi')&&strcmpi(op_req,'st')
    tempX = bi2de(ip,'left-msb');
    op = ['x' sprintf('%07d',tempX) 'x'];
elseif strcmpi(ip_type,'de')&&strcmpi(op_req,'st')
    op = ['x' sprintf('%07d',ip) 'x'];
end
end