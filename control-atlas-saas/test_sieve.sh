#!/bin/bash
API_URL="http://localhost:8000"
EMAIL="test@example.com"
PASSWORD="physics_is_law"

echo "👤 Logging in..."
TOKEN=$(curl -s -X POST "$API_URL/auth/login" -F "username=$EMAIL" -F "password=$PASSWORD" | grep -oP '(?<="access_token":")[^"]*')

echo -e "\n🧬 4. Uploading VALID minimal CIF..."
cat > valid_test.cif << 'CIF'
data_test
_entry.id test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N ALA A 1 27.340 24.430 2.614 1.00 9.67
ATOM 2 CA CA ALA A 1 28.100 25.600 3.100 1.00 9.67
ATOM 3 N N ALA B 1 35.340 30.430 8.614 1.00 9.67
ATOM 4 CA CA ALA B 1 36.100 31.600 9.100 1.00 9.67
CIF

JOB_ID=$(curl -s -X POST "$API_URL/jobs/upload" \
     -H "Authorization: Bearer $TOKEN" \
     -F "target_chain=A" \
     -F "binder_chain=B" \
     -F "file=@valid_test.cif" | grep -oP '(?<="job_id":")[^"]*')

echo "Job ID: $JOB_ID"
echo "⏳ Waiting 5s..."
sleep 5
curl -s -X GET "$API_URL/jobs/status/$JOB_ID" -H "Authorization: Bearer $TOKEN"
