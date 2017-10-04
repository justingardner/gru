function passwordless = testPasswordlessSSH( hostString )

ret = system(sprintf('ssh -o ''PreferredAuthentications=publickey'' %s "echo"',hostString));

if ret>0
    disp(sprintf('Passwordless SSH is not set up for %s: attempting to do automatic setup',hostString));
    
    system(sprintf('ssh %s; touch ~/.ssh/authorized_keys',hostString));
    system(sprintf('cat ~/.ssh/id_rsa.pub | ssh %s ''cat >> ~/.ssh/authorized_keys''',hostString));
    
    ret = system(sprintf('ssh -o ''PreferredAuthentications=publickey'' %s "echo"',hostString));
end

passwordless = ret==0;