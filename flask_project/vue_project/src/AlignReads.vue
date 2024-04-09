<script setup>
import { ref } from 'vue';
import axios from 'axios';
const file = ref(null);
const percentage = ref(0);
const task_id = ref(null);
const status = ref(null);
const ref1 = ref(null);
const ref2 = ref(null);
const cut1 = ref(null);
const cut2 = ref(null);
const PAM1 = ref("NGG");
const PAM2 = ref("NGG");

async function alignReads(event) {
    file.value = event.target.files;
    if (!file.value) {
        alert("file is not selected");
        return;
    }
    try{
        status.value = "Uploading";
        const response = await axios.putForm("/align", {
            files: file.value,
            ref1: ref1.value,
            ref2: ref2.value,
            cut1: cut1.value,
            cut2: cut2.value,
            PAM1: PAM1.value,
            PAM2: PAM2.value,
        },
        {
            onUploadProgress : progressEvent => {
                percentage.value = Math.round((progressEvent.loaded * 100) / progressEvent.total);
            }
        });
        task_id.value = response.data.split(" ")[0]
        status.value = response.data.split(" ")[1]
    } catch(error) {
        alert(error);
    }
}

async function inspect() {
    const response = await axios.get("/inspect/" + task_id.value)
    status.value = response.data
    if (status.value == 'SUCCESS') {
        const response = await axios.get("/download/" + task_id.value, {
            responseType: 'blob',
        })
        const link = document.createElement('a');
        link.href = URL.createObjectURL(new Blob([response.data]));
        link.download = response.headers['content-disposition'].split("filename=")[1];
        document.body.appendChild(link);
        link.click();
        URL.revokeObjectURL(link.href);
        document.body.removeChild(link);
    }
}
</script>

<template>
    <div>
        <input type="file" @change="alignReads" :disabled="!ref1 || !ref2 || !cut1 || !cut2"><br>
        <span>{{ percentage }}%</span><br>
        <span>{{ task_id }}</span><br>
        <span>{{ status }}</span><br>
        <button @click="inspect" :disabled="percentage < 100 || !task_id || status == 'FAILURE'">inspect</button><br>
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
