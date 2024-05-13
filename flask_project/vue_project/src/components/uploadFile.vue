<script setup>
import { ref } from 'vue';
import axios from 'axios';
const file = ref(null);
const percentage = ref(0);
const props = defineProps({
    fileType: String
});
const pageStatus = defineModel()
async function uploadFile(event) {
    file.value = event.target.files;
    try{
        const respones = await axios.putForm("/uploadFile", {
            file: file.value,
        },
        {
            onUploadProgress : progressEvent => {
                percentage.value = Math.round((progressEvent.loaded * 100) / progressEvent.total);
            }
        });
        pageStatus.value[props.fileType] = respones.data
    } catch(error) {
        alert(error);
    }
}
</script>

<template>
    <div>
        <span>{{ fileType }}</span>
        <input type="file" @change="uploadFile">
        <span>{{ percentage }}%</span>
    </div>
</template>