<script setup>
import { ref } from 'vue';
import axios from 'axios';
const file = ref(null);
const percentage = ref(0);
const ref1 = ref(null);
const ref2 = ref(null);
const cut1 = ref(null);
const cut2 = ref(null);
const PAM1 = ref("NGG");
const PAM2 = ref("NGG");

async function uploadFile(event) {
    file.value = event.target.files;
    if (!file.value) {
        alert("file is not selected");
        return;
    }
    try{
        const response = await axios.putForm("/upload", {
            files: file.value,
            ref1: ref1.value,
            ref2: ref2.value,
            cut1: cut1.value,
            cut2: cut2.value,
            PAM1: PAM1.value,
            PAM2: PAM2.value,
        },
        {
            responseType: 'blob',
            onUploadProgress : progressEvent => {
                percentage.value = Math.round((progressEvent.loaded * 100) / progressEvent.total);
            }
        });
        const link = document.createElement('a');
        link.href = URL.createObjectURL(new Blob([response.data]));
        link.download = response.headers['content-disposition'].split("filename=")[1];
        document.body.appendChild(link);
        link.click();
        URL.revokeObjectURL(link.href);
        document.body.removeChild(link);
    } catch(error) {
        alert(error);
    }
}
</script>

<template>
    <div>
        <input type="file" @change="uploadFile" :disabled="!ref1 || !ref2 || !cut1 || !cut2"><br>
        <span>"{{ percentage }}%"</span><br>
        <textarea v-model="ref1" placeholder="reference1"></textarea><br>
        <textarea v-model="ref2" placeholder="reference2"></textarea><br>
        <input v-model="cut1" placeholder="cleavage1">
        <input v-model="cut2" placeholder="cleavage2"><br>
        <select v-model="PAM1" :required>
            <option>NGG</option>
            <option>CCN</option>
        </select>
        <select v-model="PAM2" :required>
            <option>NGG</option>
            <option>CCN</option>
        </select>
    </div>
</template>