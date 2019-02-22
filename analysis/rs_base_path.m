% Set up paths for windows vs unix
host_name = char(java.net.InetAddress.getLocalHost.getHostName);
if strcmp(computer, 'PCWIN64')
    base_dir = 'Z:/gb/';
elseif strcmp(host_name, 'colles-d164179')
    base_dir = '/home/local/ADF/brookshg/rds_share/gb/';
else
    base_dir = '/rds/projects/2017/jenseno-02/gb/';
end